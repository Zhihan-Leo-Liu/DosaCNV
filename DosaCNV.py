import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input, Dense, Lambda, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras import regularizers
import time
import argparse
import random
import os

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest = 'subparser_name')

  parser_a = subparsers.add_parser("train")
  parser_a.add_argument("-train_data", dest="train_index", required=True,
                        help="input training data")
  parser_a.add_argument("-val_data", dest="val_index", required=True,
                        help="input validation data")
  parser_a.add_argument("-annotation", dest="annotation", required=True,
                        help="gene level annotation")
  parser_a.add_argument("-max_n_gene", dest="max_n_gene", default = 100, type = int,
                        help="maximum number of genes within a deletion")
  parser_a.add_argument("-output_model_name", dest="output_model_name", required=True,
                        help="set name for saved model/models")
  parser_a.add_argument("-lr", dest="learning_rate", default = 3e-4, type = float,
                        help="model learning rate")
  parser_a.add_argument("-seed", dest="seed",default = int(time.time()), type = int,
                        help="set seed for training")
  parser_a.add_argument("-save_gene_model", dest="save_gene_model",action='store_true',
                        help="save gene level model (DosaCNV-HI)")

  parser_b = subparsers.add_parser("predict-deletion")
  parser_b.add_argument("-input_data", dest="input_index", required=True,
                        help="input deletion data")
  parser_b.add_argument("-annotation", dest="annotation", default="./annotation/anno_for_deletion.txt",
                        help="gene level annotation")
  parser_b.add_argument("-saved_model_name", dest="saved_model_name",default="pretrained_deletion",
                        help="load saved model for prediction")
  parser_b.add_argument("-output_name", dest="output_name", required = True,
                        help="output file name")
  parser_b.add_argument("-seed", dest="seed",default = int(time.time()), type = int,
                        help="set seed")
  parser_b.add_argument("-max_n_gene", dest="max_n_gene", default = 100, type = int,
                        help="maximum number of genes within a deletion")

  parser_c = subparsers.add_parser("predict-gene")
  parser_c.add_argument("-input_data", dest="input_index", required=True,
                        help="input deletion data")
  parser_c.add_argument("-saved_model_name", dest="saved_model_name", default="pretrained_gene",
                        help="load saved model for prediction")
  parser_c.add_argument("-output_name", dest="output_name", required = True,
                        help="output file name")
  parser_c.add_argument("-seed", dest="seed",default = int(time.time()), type = int,
                        help="set seed")

  parser_d = subparsers.add_parser("explain-gene")
  parser_d.add_argument("-explain_set", dest="explain_set", required=True,
                        help="explaining set")
  parser_d.add_argument("-background_set", dest="background_set", default="./data/1000_hs_gene.txt",
                        help="background set")
  parser_d.add_argument("-saved_model_name", dest="saved_model_name", default="pretrained_gene", 
                        help="load saved model to initialize explainer")
  parser_d.add_argument("-output_name", dest="output_name", required = True,
                        help="output file name")
  parser_d.add_argument("-seed", dest="seed",default = int(time.time()), type = int,
                        help="set seed")

  args = parser.parse_args()

subparser_name = vars(args)['subparser_name']

seed_value= args.seed
os.environ['PYTHONHASHSEED']=str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)
tf.random.set_seed(seed_value)
from keras import backend as K
session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
tf.compat.v1.keras.backend.set_session(sess)

def annot(data,annotation):
  combined = data.merge(annotation, how='left', on='gene_id')
  if 'label' in combined.columns:
    x = combined.drop(['gene_id','label','type'],axis =1)
    y = combined.drop_duplicates(subset = ['variant_id','label','type'])
  else:
    x = combined.drop(['gene_id','type'],axis =1)
    y = combined.drop_duplicates(subset = ['variant_id','type'])
  return (x,y)

def padding(data):
  padded_inputs = []
  num_col= len(data.columns)
  num_feat = num_col-1
  for id, df_id in data.groupby('variant_id',sort=False):
    df =df_id.to_numpy(copy=True)
    df1 = df.reshape(1,-1,num_col)
    df2 = df1[:,:,1:]
    tmp_padded_inputs = pad_sequences(df2, padding="post", maxlen=args.max_n_gene, dtype='float32')
    padded_inputs.append(tmp_padded_inputs)
  padded_inputs = np.vstack(padded_inputs)
  return padded_inputs

if not os.path.exists('./training_history/'):
  os.makedirs('./training_history/')
if not os.path.exists('./variant_score/'):
  os.makedirs('./variant_score/')
if not os.path.exists('./gene_score/'):
  os.makedirs('./gene_score/')
if not os.path.exists('./shap_value/'):
  os.makedirs('./shap_value/')

def training(x, y, val_x, val_y,
             model_name,
             save_gene_model = False,
             max_n_gene=100,
             units_1=128,
             units_2=64,
             l2_reg=0.02,
             dropout_rate=0.5,
             learning_rate=3e-4,
             batch_size=32,
             epochs=1000,
             patience=15):

  num_feat = len(x.columns)-1
  dtype = 'float32'
  val_x_padded = padding(val_x)
  train_x_padded = padding(x)

  #DosaCNV main model
  feature_input = Input(shape=(max_n_gene,num_feat),name="input_layer",dtype=dtype)

  def generate_mask(x):
    mask = tf.math.not_equal(tf.reduce_sum(x,axis=2),0)
    return mask

  gene_mask = Lambda(generate_mask,name='gene_mask')(feature_input)

  hidden1 = Dense(units=units_1, name='hidden1', activation='relu',kernel_regularizer=regularizers.l2(l2 = l2_reg),
                 kernel_initializer='he_normal')(feature_input)

  hidden2 =Dense(units=units_2, name='hidden2', activation='relu',kernel_regularizer=regularizers.l2(l2 = l2_reg),
                 kernel_initializer='he_normal')(hidden1)

  drop = Dropout(dropout_rate ,name='drop')(hidden2)

  gene_HI_risk = Dense(units=1, name='gene_HI_risk',activation='sigmoid')(drop)

  HI_gene_predictor = Model(feature_input, gene_HI_risk, name = 'DosaCNV-HI')

  def noisy_or_function(x):
    gene_risk = tf.squeeze(x[0])
    mask = tf.cast(x[1], dtype)
    total_risk = 1 - tf.math.reduce_prod(1 - gene_risk*mask, axis=1)
    return total_risk

  total_risk = Lambda(noisy_or_function, name='total_risk')([gene_HI_risk, gene_mask])

  model = Model(inputs=feature_input, outputs= total_risk, name="DosaCNV")

  model.compile(optimizer = Adam(learning_rate=learning_rate),
                loss= tf.keras.losses.BinaryCrossentropy(),metrics=['binary_accuracy','AUC'])

  callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                            restore_best_weights=True,
                                            patience=patience)

  history = model.fit(train_x_padded,y['label'].values, epochs=epochs,
                      batch_size=batch_size, verbose = 1,
                      validation_data = (val_x_padded,val_y['label'].values),
                      callbacks=[callback])

  pd.DataFrame(history.history).to_csv(r'' + './training_history/' + model_name
                                       + '_training_history.tsv',index=False, sep='\t',header = True)


  model.save(r'' + './saved_model/' + model_name +'_deletion')  #DosaCNV: variant-level model


  #extract all layers from the main model except for the noisy-or pooling layer
  config = HI_gene_predictor.get_config()     #get the current model configuration

  input_shape = config['layers'][0]['config']['batch_input_shape']       #modify the configuration (batch_size,max_n_gene,num_feat) to

  new_input_shape = (input_shape[0], input_shape[-1])                    #(batch_size,num_feat)

  config['layers'][0]['config']['batch_input_shape'] = new_input_shape

  hi_gene_score = tf.keras.Model().from_config(config)   #create a new model from the modified configuration

  for i in range(len(HI_gene_predictor.layers)):
    hi_gene_score.layers[i].set_weights(HI_gene_predictor.layers[i].get_weights())

  if save_gene_model == True:
    hi_gene_score.save(r'' + './saved_model/' + model_name + '_gene')  # DosaCNV-HI: gene-level model

  return (model,hi_gene_score)

# training and saving variant and/or gene level models
if subparser_name == 'train':
  train_index = pd.read_csv(args.train_index,sep=',')
  val_index = pd.read_csv(args.val_index,sep=',')
  annotation =  pd.read_csv(args.annotation, sep=',')

  train_x,train_y = annot(train_index,annotation)
  val_x,val_y = annot(val_index,annotation)

  model = training(train_x, train_y, val_x, val_y,
                   model_name = args.output_model_name,
                   save_gene_model = args.save_gene_model,
                   max_n_gene = args.max_n_gene,
                   learning_rate = args.learning_rate)

# variant-level prediction
elif subparser_name == 'predict-deletion':

  del_index = pd.read_csv(args.input_index,sep=',')
  annotation = pd.read_csv(args.annotation, sep=',')
  x,y = annot(del_index,annotation)

  del_model = tf.keras.models.load_model(r'' + './saved_model/' + args.saved_model_name)
  del_score = del_model.predict(padding(x))

  del_id = x.drop_duplicates(subset = ['variant_id'])['variant_id']
  del_id = del_id.reset_index()
  del_id = del_id.drop(['index'],axis =1)
  prediction = pd.concat([del_id,pd.DataFrame(del_score,columns = ['DosaCNV_score'])],axis = 1)

  pd.DataFrame(prediction).to_csv(r''+'./variant_score/'+ args.output_name +
                                  '_score.tsv',index=False, sep='\t',header = True)

# gene-level prediction
elif subparser_name == 'predict-gene':
  gene_input = pd.read_csv(args.input_index, sep=',')
  gene_feature = gene_input.loc[:,gene_input.columns!= 'gene_id'].to_numpy()
  gene_id = gene_input['gene_id']

  gene_model = tf.keras.models.load_model(r'' + './saved_model/' + args.saved_model_name)

  gene_score = gene_model.predict(gene_feature)
  gene_score = pd.concat([gene_id,pd.DataFrame(gene_score,columns=['DosaCNV_HI_score'])],axis=1)
  gene_score.to_csv(r'' + './gene_score/' + args.output_name + '_score.tsv',index=False, sep='\t',header = True)

# shap values for each gene across all features
elif subparser_name == 'explain-gene':
  try:
    import shap
  except ImportError:
    print("The 'shap' package is not installed. Please install it by running 'pip install shap' in your terminal.")
    import sys
    sys.exit(1)

  explain_set = pd.read_csv(args.explain_set, sep=',')
  explain_set_feature = explain_set.loc[:,explain_set.columns!= 'gene_id'].to_numpy()
  gene_feat_names = explain_set.loc[:,explain_set.columns!= 'gene_id'].columns
  explain_gene_id = explain_set['gene_id']

  background_set = pd.read_csv(args.background_set, sep=',')
  background_set_feature = background_set.loc[:,background_set.columns!= 'gene_id'].to_numpy()
  background_gene_id = background_set['gene_id']

  gene_model = tf.keras.models.load_model(r'' + './saved_model/' + args.saved_model_name)

  explainer = shap.DeepExplainer(gene_model,background_set_feature)
  shap_values = explainer.shap_values(explain_set_feature)
  shap_values = np.array(shap_values).reshape((-1,len(explain_set.columns)-1))
  shap_values_df = pd.concat([explain_gene_id,
                            pd.DataFrame(shap_values,columns=gene_feat_names)],axis=1)

  shap_values_df.to_csv(r'' + './shap_value/' + args.output_name + '_shap.csv',index=False,sep=',',header=True)