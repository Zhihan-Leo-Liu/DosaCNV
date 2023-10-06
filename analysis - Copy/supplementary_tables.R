library(data.table)
library(MatchIt)
library(pROC)
library(ggplot2)
library(reshape)
library(scales)
library(cowplot)
library(dplyr)

######################
#Supplementary table 1
######################

our_nstd102_scores <- fread('./our_scores/nstd102_score.tsv')
colnames(our_nstd102_scores)[2] <- 'our_score'
meta_nstd102_scores <- fread('./meta_scores/meta_nstd102_scores.csv')

meta_nstd102_scores$var_length <- meta_nstd102_scores$end_hg19 - meta_nstd102_scores$start_hg19

combined_variant <- merge(our_nstd102_scores,meta_nstd102_scores,by = 'variant_id',all.x = T, all.y = F)
combined_variant <- combined_variant[complete.cases(combined_variant),]
combined_variant <- combined_variant[,c(2,11,12,13,14,15,17,19,18)]

seed = 1

var_comp <- function(scores, quiet = F) {
  # Ensure scores is a dataframe
  scores <- as.data.frame(scores)
  
  # Display basic information unless quiet is TRUE
  if (!quiet) {
    cat("samples:", nrow(scores), "mean:", mean(scores$label), "\n")
  }
  
  # Identify number of methods and their names
  num_methods <- ncol(scores) - 1
  methods_name <- colnames(scores)[1:num_methods]
  
  # Compute ROC object
  roc_obj <- roc(scores$label ~ ., data = scores, quiet = TRUE)
  
  # Compute AUC values and round to 4 decimal places
  auc_value <- round(sapply(roc_obj, function(x) x$auc), 4)
  
  # Compute SE for AUC
  se_value <- round(sapply(roc_obj, function(x) as.vector(ci(x))), 4)   ###
  se_value <- (se_value[3,] - se_value[1,]) / (1.96 * 2)                ###
  se_value <- round(se_value, 3)                                        ###
  
  # Rank methods based on AUC
  methods_rank <- names(auc_value[order(-auc_value)]) 
  
  # Compute p-value for 1st and 2nd ranked methods
  pvalue <- roc.test(roc_obj[[methods_rank[1]]],roc_obj[[methods_rank[2]]])$p.value  ###

  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, pvalue), nrow = 1))  ###
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results unless quiet is TRUE
  if (!quiet) {
    print(list_result)
  }
  
  return(list(roc_obj, list_result, se_value))  ###
}

var_comp_match_both <- function(scores) {
  # Convert scores to a dataframe
  scores <- as.data.frame(scores)
  
  # Display number of samples before matching
  cat("samples before matching:", nrow(scores), "\n")
  
  # Matching process
  set.seed(seed)
  mout <- matchit(
    formula = label ~ var_length,
    data = scores, 
    exact = ~num_coding_gene,
    method = 'nearest',
    caliper = 1,
    discard = 'both',
    m.order = 'random'
  )
  scores <- match.data(mout)
  
  # Remove unnecessary columns after matching
  scores <- subset(scores, select = -c((ncol(scores)-2):ncol(scores)))
  
  # Display number of samples after matching
  cat("samples after matching:", nrow(scores), "mean:", mean(scores$label), "\n")
  
  # Identify number of methods and their names
  num_methods <- ncol(scores) - 1
  methods_name <- colnames(scores)[1:num_methods]
  
  # Compute ROC object
  roc_obj <- roc(scores$label ~ ., data = scores, quiet = TRUE)
  
  # Compute AUC values and round to 4 decimal places
  auc_value <- round(sapply(roc_obj, function(x) x$auc), 4)
  
  # Compute SE for AUC
  se_value <- round(sapply(roc_obj, function(x) as.vector(ci(x))), 4)   ###
  se_value <- (se_value[3,] - se_value[1,]) / (1.96 * 2)                ###
  se_value <- round(se_value, 3)                                        ###
  
  # Rank methods based on AUC (this isn't used later, but kept for consistency)
  methods_rank <- names(auc_value[order(-auc_value)])
  
  # Compute p-value for our_score vs. tada
  pvalue <- roc.test(roc_obj[[methods_rank[1]]],roc_obj[[methods_rank[2]]])$p.value  ###
  
  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, pvalue), nrow = 1))  ###
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results
  print(list_result)
  
  return(list(roc_obj, list_result, se_value))  ###
}

var_comp_match_cdg <- function(scores) {
  # Convert scores to a dataframe
  scores <- as.data.frame(scores)
  
  # Display number of samples before matching
  cat("samples before matching:", nrow(scores), "\n")
  
  # Matching process
  set.seed(seed)
  mout <- matchit(
    formula = label ~ num_coding_gene,
    data = scores,
    method = 'nearest',
    caliper = 0,
    discard = 'both',
    m.order = 'random'
  )
  scores <- match.data(mout)
  
  # Remove unnecessary columns after matching
  scores <- subset(scores, select = -c((ncol(scores)-2):ncol(scores)))
  
  # Display number of samples after matching
  cat("samples after matching:", nrow(scores), "mean:", mean(scores$label), "\n")
  
  # Identify number of methods and their names
  num_methods <- ncol(scores) - 1
  methods_name <- colnames(scores)[1:num_methods]
  
  # Compute ROC object
  roc_obj <- roc(scores$label ~ ., data = scores, quiet = TRUE)
  
  # Compute AUC values and round to 4 decimal places
  auc_value <- round(sapply(roc_obj, function(x) x$auc), 4)
  
  # Compute SE for AUC
  se_value <- round(sapply(roc_obj, function(x) as.vector(ci(x))), 4)   ###
  se_value <- (se_value[3,] - se_value[1,]) / (1.96 * 2)                ###
  se_value <- round(se_value, 3)                                        ###
  
  # Rank methods based on AUC (this isn't used later, but kept for consistency)
  methods_rank <- names(auc_value[order(-auc_value)])
  
  # Compute p-value for our_score vs. tada
  pvalue <- roc.test(roc_obj[[methods_rank[1]]],roc_obj[[methods_rank[2]]])$p.value  ###
  
  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, pvalue), nrow = 1))  ###
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results
  print(list_result)
  
  return(list(roc_obj, list_result, se_value))
}

var_comp_match_length <- function(scores) {
  # Convert scores to a dataframe
  scores <- as.data.frame(scores)
  
  # Display number of samples before matching
  cat("samples before matching:", nrow(scores), "\n")
  
  # Matching process
  set.seed(seed)
  mout <- matchit(
    formula = label ~ var_length,
    data = scores,
    method = 'nearest',
    caliper = 1,
    discard = 'both',
    m.order = 'random'
  )
  scores <- match.data(mout)
  
  # Remove unnecessary columns after matching
  scores <- subset(scores, select = -c((ncol(scores)-2):ncol(scores)))
  
  # Display number of samples after matching
  cat("samples after matching:", nrow(scores), "mean:", mean(scores$label), "\n")
  
  # Identify number of methods and their names
  num_methods <- ncol(scores) - 1
  methods_name <- colnames(scores)[1:num_methods]
  
  # Compute ROC object
  roc_obj <- roc(scores$label ~ ., data = scores, quiet = TRUE)
  
  # Compute AUC values and round to 4 decimal places
  auc_value <- round(sapply(roc_obj, function(x) x$auc), 4)
  
  # Compute SE for AUC
  se_value <- round(sapply(roc_obj, function(x) as.vector(ci(x))), 4)   ###
  se_value <- (se_value[3,] - se_value[1,]) / (1.96 * 2)                ###
  se_value <- round(se_value, 3)                                        ###
  
  # Rank methods based on AUC (this isn't used later, but kept for consistency)
  methods_rank <- names(auc_value[order(-auc_value)])
  
  # Compute p-value for our_score vs. tada
  pvalue <- roc.test(roc_obj[[methods_rank[1]]],roc_obj[[methods_rank[2]]])$p.value  ###
  
  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, pvalue), nrow = 1))  ###
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results
  print(list_result)
  
  return(list(roc_obj, list_result, se_value))  ###
}

size_group_comp <- function(scores, matchornot) {
  
  # Convert scores to a dataframe
  scores <- as.data.frame(scores)
  set.seed(2)
  
  # Match scores based on the 'matchornot' parameter
  if (matchornot == 'cdg') {
    mout <- matchit(formula = label ~ num_coding_gene, data = scores,
                    caliper = 0, method = 'nearest',
                    discard = 'both', m.order = 'random')
  } else if (matchornot == 'length') {
    mout <- matchit(formula = label ~ var_length, data = scores,
                    method = 'nearest', caliper = 1,
                    discard = 'both', m.order = 'random')
  } else if (matchornot == 'both') {
    mout <- matchit(formula = label ~ var_length, data = scores, exact = ~num_coding_gene,
                    method = 'nearest', caliper = 1,
                    discard = 'both', m.order = 'random')
  }
  
  if (matchornot != 'unmatch') {
    scores <- match.data(mout)
    scores <- subset(scores, select = -c((ncol(scores) - 2):ncol(scores)))
  }
  
  
  # Define size groups based on var_length
  q1 <- scores[scores$var_length <= 50000, ]
  q2 <- scores[scores$var_length > 50000 & scores$var_length <= 100000, ]
  q3 <- scores[scores$var_length > 100000 & scores$var_length <= 1000000, ]
  q4 <- scores[scores$var_length > 1000000, ]
  
  # Extract AUC values for each size group
  q1_auc <- var_comp(q1, quiet = TRUE)[[2]]
  q2_auc <- var_comp(q2, quiet = TRUE)[[2]]
  q3_auc <- var_comp(q3, quiet = TRUE)[[2]]
  q4_auc <- var_comp(q4, quiet = TRUE)[[2]]
  
  out <- rbind.data.frame(q1_auc, q2_auc, q3_auc, q4_auc)
  out <- out[, -10]
  
  # Helper function to extract AUC information
  extract_auc_data <- function(data) {
    t(as.data.frame(sapply(var_comp(data, quiet = TRUE)[[1]], function(x) ci(x)))[, -10])
  }
  
  q1_auc_full <- extract_auc_data(q1)
  q2_auc_full <- extract_auc_data(q2)
  q3_auc_full <- extract_auc_data(q3)
  q4_auc_full <- extract_auc_data(q4)
  
  q1_se <- format(round((q1_auc_full[,3]-q1_auc_full[,1])/3.92,3),nsmall=3)  
  q2_se <- format(round((q2_auc_full[,3]-q2_auc_full[,1])/3.92,3),nsmall=3)
  q3_se <- format(round((q3_auc_full[,3]-q3_auc_full[,1])/3.92,3),nsmall=3)
  q4_se <- format(round((q4_auc_full[,3]-q4_auc_full[,1])/3.92,3),nsmall=3)
  se_value <- rbind.data.frame(q1_se,q2_se,q3_se,q4_se)
  colnames(se_value) <- NULL
  
  out_plot <- rbind.data.frame(q1_auc_full, q2_auc_full, q3_auc_full, q4_auc_full)
  out_plot$q <- rep(c('<50 ', '50-100 ', '100-1000 ', '>1000 '), each = nrow(out_plot) / 4)
  out_plot$method <- rep(colnames(scores)[1:8], 4)
  colnames(out_plot)[1:3] <- c('lower', 'average', 'upper')
  out_plot$upper <- out_plot$average + (out_plot$upper - out_plot$average) / 1.96 
  out_plot$lower <- out_plot$average - (out_plot$average - out_plot$lower) / 1.96
  
  num_q <- c(nrow(q1), nrow(q2), nrow(q3), nrow(q4))
  
  print(out)
  
  return(list(out, out_plot, num_q, se_value))
}

compute_var_stats <- function(function_to_call, variant_data) {
  result <- function_to_call(variant_data)
  
  auc_vals <- round(result[[2]][, 1:8], 3)
  se_vals <- result[[3]]
  formatted_auc_vals <- paste(format(auc_vals, nsmall=3), "±", format(se_vals, nsmall=3), sep="")
  p_value <- format(result[[2]]$pvalue, nsmall=1)
  
  c(formatted_auc_vals, p_value)
}

# Compute stats for different methods
stats_default <- compute_var_stats(var_comp, combined_variant)
stats_match_cdg <- compute_var_stats(var_comp_match_cdg, combined_variant)
stats_match_length <- compute_var_stats(var_comp_match_length, combined_variant)
stats_match_both <- compute_var_stats(var_comp_match_both, combined_variant)

method_rank <- c('DosaCNV', 'TADA', 'StrVCTVRE', 'HIS-LOD', 'X-CNV', 'CADD-SV', 'num affected gene', 'deletion length', 'p-value')
method_names <- c('DosaCNV', 'CADD-SV', 'TADA', 'X-CNV', 'StrVCTVRE', 'HIS-LOD', 'num affected gene', 'deletion length', 'p-value')
combined_stats <- rbind.data.frame(method_names, stats_default, stats_match_cdg, stats_match_length, stats_match_both)

table_s1 <- as.data.frame(t(combined_stats))
row.names(table_s1) <- NULL
table_s1 <- table_s1[match(method_rank,table_s1$V1),]
colnames(table_s1) <- c('','Unmatched','Matched by number of affected genes','Matched by deletion length','Matched by both')
write.table(table_s1,'table_s1.csv',quote = F,row.names = F,sep=',')

######################
#Supplementary table 2 & 3
######################

# Function to process and format the data for the specified condition
process_var_data <- function(variant_data, condition, method_order) {
  data_list <- size_group_comp(variant_data, condition)
  
  formatted_data <- method_names
  for (i in 1:4) {
    tmp <- c(paste(format(round(data_list[[1]][i, ][1:8], 3), nsmall = 3), 
                   '±', data_list[[4]][i, ], sep = ''), data_list[[1]][i, 9])
    formatted_data <- rbind.data.frame(formatted_data, tmp)
  }
  
  formatted_data <- as.data.frame(t(formatted_data))
  formatted_data <- formatted_data[match(method_order, formatted_data$V1),]
  row.names(formatted_data) <- NULL
  
  if (condition == "unmatch") {
    formatted_data$matching_condition <- 'Unmatched'
  } else if (condition == "both") {
    formatted_data$matching_condition <- 'Both'
  } else if (condition == "cdg") {
    formatted_data$matching_condition <- 'Number of affected genes'
  } else {
    formatted_data$matching_condition <- 'Deletion length'
  }
  
  return(formatted_data)
}

# Function to interleave two data frames
interleave_var_data <- function(df1, df2) {
  interleaved <- rbind(df1[1, ], df2[1, ])
  for (i in 2:nrow(df1)) {
    interleaved <- rbind(interleaved, df1[i, ], df2[i, ])
  }
  
  # Reorder and rename columns
  interleaved <- interleaved[, c(1, 6, 2, 3, 4, 5)]
  colnames(interleaved) <- c('', 'Matching condition', 'Small (<50 kb)', 
                             'Medium (50-100 kb)', 'Large (100-1000 kb)', 
                             'Extra-large (>1000 kb)')
  
  return(interleaved)
}

# Process data for table_s2
unmatched_data <- process_var_data(combined_variant, 'unmatch', method_rank)
both_data <- process_var_data(combined_variant, 'both', method_rank)
table_s2 <- interleave_var_data(unmatched_data, both_data)

# Process data for table_s3
cdg_data <- process_var_data(combined_variant, 'cdg', method_rank)
length_data <- process_var_data(combined_variant, 'length', method_rank)
table_s3 <- interleave_var_data(cdg_data, length_data)

write.table(table_s2,'table_s2.csv',quote = F,row.names = F,sep=',')
write.table(table_s3,'table_s3.csv',quote = F,row.names = F,sep=',')

######################
#Supplementary table 4
######################

our_nstd100_score <- fread('./our_scores/nstd100_score.tsv')

meta_nstd100_score <- fread('./meta_scores/meta_nstd100_scores.csv')

our_nstd173_score <- fread('./our_scores/nstd173_score.tsv')

meta_nstd173_score <- fread('./meta_scores/meta_nstd173_scores.csv')

colnames(our_nstd100_score)[2] <- 'our_score'

colnames(our_nstd173_score)[2] <- 'our_score'

case_control_enrich <- function(meta, our, study, section = "top") {
  meta <- as.data.frame(meta)
  our <- as.data.frame(our)
  # Determine the number of tops based on the study type
  if (study == 'nstd100') {
    num_tops <- 1454
  } else if (study == 'nstd173') {
    num_tops <- 134
  } else {
    stop('Invalid study type. Please specify a valid study.') 
  }
  
  # Merge datasets and select necessary columns
  combined <- merge(meta, our, by = 'variant_id')
  combined <- combined[, c(9:15, 18, 16:17)]
  
  # Filter rows based on num_coding_gene
  combined <- combined[combined$num_coding_gene <= 100, ]
  
  method_names <- colnames(combined)[2:9]
  print(method_names)
  
  case_counts <- c()
  control_counts <- c()
  
  # Calculate case and control counts for each method
  for (i in 2:10) {
    data_ordered <- combined[order(combined[, i], decreasing = TRUE), ]
    
    # Select the top or bottom data based on the section parameter
    if (section == "top") {
      selected_data <- data_ordered[1:num_tops, ]
    } else if (section == "bottom") {
      selected_data <- data_ordered[num_tops:nrow(data_ordered), ]
    } else {
      stop("Invalid section type. Please choose either 'top' or 'bottom'.")
    }
    
    case_counts <- c(case_counts, sum(selected_data$status == 1))
    control_counts <- c(control_counts, sum(selected_data$status == 0))
  }
  
  # Convert counts to data frames
  case_counts_df <- as.data.frame(t(matrix(case_counts, nrow = 1)))
  control_counts_df <- as.data.frame(t(matrix(control_counts, nrow = 1)))
  
  # Combine counts into a single dataframe and set method names
  combined_counts <- data.frame(
    case_counts = case_counts_df,
    control_counts = control_counts_df,
    method = c('CADD-SV', 'StrVCTVRE', 'AnnotSV', 'TADA', 'X-CNV', 'HIS-LOD', 'DosaCNV', 'num affected gene', 'deletion length'),
    study = study
  )
  
  # Name columns based on the section parameter
  if (section == "top") {
    colnames(combined_counts)[1:2] <- c('case_counts_top', 'control_counts_top')
  } else {
    colnames(combined_counts)[1:2] <- c('case_counts_bot', 'control_counts_bot')
  }
  
  # Order factors for the method
  combined_counts$method <- factor(combined_counts$method, levels = combined_counts$method[order(combined_counts$case_counts)])
  
  return(combined_counts)
}

calculate_or_and_se <- function(data) {
  data$or <- log((data$case_counts_top / data$control_counts_top) / 
                   (data$case_counts_bot / data$control_counts_bot))
  data$or_se <- sqrt(1 / data$case_counts_top + 1 / data$control_counts_top + 
                       1 / data$case_counts_bot + 1 / data$control_counts_bot)
  return(data)
}

nstd100_data <- case_control_enrich(meta_nstd100_score, our_nstd100_score, 'nstd100', 'top')
nstd173_data <- case_control_enrich(meta_nstd173_score, our_nstd173_score, 'nstd173', 'top')

nstd100_data_bottom <- case_control_enrich(meta_nstd100_score, our_nstd100_score, 'nstd100', 'bottom')
nstd173_data_bottom <- case_control_enrich(meta_nstd173_score, our_nstd173_score, 'nstd173', 'bottom')

# Merging top and bottom data
nstd100_data <- merge(nstd100_data, nstd100_data_bottom, by=c('method', 'study'))
nstd173_data <- merge(nstd173_data, nstd173_data_bottom, by=c('method', 'study'))

# Calculating OR and SE
nstd100_data <- calculate_or_and_se(nstd100_data)
nstd173_data <- calculate_or_and_se(nstd173_data)

# Calculating upper bound for plotting
nstd100_data$upper <- exp(nstd100_data$or + nstd100_data$or_se)
nstd173_data$upper <- exp(nstd173_data$or + nstd173_data$or_se)

# P values for difference between DosaCNV and HIS-LOD
nstd100_zscore <- (4.627555-3.436920)/sqrt(0.0951425^2 + 0.1577527^2)
nstd100_pvalue <- 2*(1-pnorm(nstd100_zscore))     # 1.0263e-10 
# P values for difference between DosaCNV and TADA
nstd173_zscore <- (0.751717890-0.583186045)/sqrt(0.2097265^2 + 0.2011201^2)
nstd173_pvalue <- 2*(1-pnorm(nstd173_zscore))     # 0.56192 

nstd100_data <- nstd100_data[-1,]  #remove annotsv
nstd173_data <- nstd173_data[-1,]  #remove annotsv

nstd100_data$true_or <- exp(nstd100_data$or)
nstd173_data$true_or <- exp(nstd173_data$or)

table_s4 <- cbind.data.frame(nstd100_data[,c('method','true_or','or','or_se')],nstd173_data[,c('true_or','or','or_se')])
table_s4 <- cbind.data.frame(table_s4$method,as.data.frame(apply(table_s4[,-1],2,function(p) format(round(p,3),nsmall=3))))
table_s4$or <- paste(table_s4$or,'±',table_s4$or_se,sep = '')
table_s4$or.1 <- paste(table_s4$or.1,'±',table_s4$or_se.1,sep = '')
table_s4 <- table_s4[,c(1,2,3,5,6)]
colnames(table_s4) <- c('method','OR','log(OR)','OR','log(OR)')
table_s4 <- table_s4[match(method_rank,table_s4$method),][1:8,]

write.table(table_s4,'table_s4.csv',quote = F,row.names = F,sep=',')

######################
#Supplementary table 5 & 6
######################

meta <- fread('./meta_scores/meta_gene_scores.csv') 
meta[meta=='missing'] <- NA 
meta <- as.data.frame(meta)
meta[,3:ncol(meta)] <- as.data.frame(apply(meta[,3:ncol(meta)],2,as.numeric))
meta <- meta[!duplicated(meta),]

our <- as.data.frame(fread('./our_scores/precomputed_gene_score.tsv'))
our <- our[,c(1,2)]
meta1 <- meta[complete.cases(meta),]
meta1 <- merge(meta1,our,by ='gene_id',all.x=T)

# number of genes within each gene list
colMeans(meta[,5:14])*nrow(meta)
colMeans(meta1[,5:14])*nrow(meta1)
colMeans(meta1[meta1$train==0,][,5:14])*nrow(meta1[meta1$train==0,])

# meta AUROC comparison
meta_auroc <- function(data){
  list_results <- c()
  lowers <- c()
  uppers <- c()
  
  for (i in c('hi','ad','mhko','ce','asd','eyeg2p','sking2p','ar','ddg2p')){
    list_name = i
    pos <- data[data[list_name]==1,]
    neg <- data[data$hs==1,]
    pos_neg <- rbind.data.frame(pos,neg)
    set.seed(1)
    mout <- matchit(formula = formula(paste(list_name,'~cds_length',collapse = '')),data=pos_neg,
                    caliper = 1,method = 'nearest',
                    discard = 'both',m.order = 'random')
    pos_neg <- match.data(mout)
    pos_neg <- subset(pos_neg,select = -c((ncol(pos_neg)-2) : ncol(pos_neg)))
    neg <- pos_neg[pos_neg[list_name]==0,]
    
    combined <- rbind.data.frame(pos,neg)
    print(paste(list_name,nrow(pos),nrow(combined)))
    combined <- cbind.data.frame(combined[list_name],combined[,17:ncol(combined)])
    roc_object <- roc(combined[,1]~.,data=combined,quiet=T)[-1]
    auc_value <- lapply(roc_object,function(x) as.vector(ci(x)))
    lower <- unlist(lapply(auc_value, function(l) l[[1]]))
    average <- unlist(lapply(auc_value, function(l) l[[2]]))
    upper <- unlist(lapply(auc_value, function(l) l[[3]]))
    
    lower <- average - (average-lower)/1.96
    upper <- (upper-average)/1.96
    
    lower <- round(lower,3)
    average <- round(average,3)
    upper <-  round(upper,3)
    
    
    methods_rank <- names(average[order(-average)])
    pvalue <- roc.test(roc_object[[methods_rank[1]]],roc_object[[methods_rank[2]]])$p.value
    list_result <- c(average,methods_rank[1:2],round(pvalue,4))
    list_results <- rbind.data.frame(list_results,list_result)
    lowers <- rbind.data.frame(lowers,lower)
    uppers <- rbind.data.frame(uppers,upper)
  }
  method_names <- c('LOEUF','Mis. OEUF','CDS length','HIS','IS','GHIS','RVIS','HIPred','VIRLoF','Episcore','EDS','pHaplo','DeepLOF',
                    'DosaCNV-HI')
  colnames(list_results) <- c(method_names,'1st','2nd','pvalue')
  row.names(list_results) <- c('hi','ad','mhko','ce','asd','eyeg2p','sking2p','ar','ddg2p')
  colnames(lowers) <- method_names
  row.names(lowers) <- c('hi','ad','mhko','ce','asd','eyeg2p','sking2p','ar','ddg2p')
  colnames(uppers) <- method_names
  row.names(uppers) <- c('hi','ad','mhko','ce','asd','eyeg2p','sking2p','ar','ddg2p')
  #print(lowers)
  print(list_results)
  #print(uppers)
  return(list(list_results,lowers,uppers))  
}

auroc_all <- meta_auroc(meta1)
auroc_test <- meta_auroc(meta1[meta1$train==0,])

table_s5_data <- auroc_all[[1]][c(1,3,5,9,6,7),]
table_s5_data1 <- auroc_all[[1]][c(1,3,5,9,6,7),][,c(1:14)]
table_s5_data1 <- data.frame(apply(table_s5_data1, 2, function(x) as.numeric(as.character(x))))
table_s5_data2 <- auroc_all[[3]][c(1,3,5,9,6,7),]
table_s5 <- c(colnames(table_s5_data2),'p-value')
for(i in 1:6){
  tmp <- c(paste(format(table_s5_data1[i,],nsmall=3),'±',format(table_s5_data2[i,],nsmall=3),sep = ''),table_s5_data[i,17])
  table_s5 <- rbind(table_s5,tmp)  
}
table_s5 <- as.data.frame(t(table_s5))
table_s5 <- table_s5[c(14,10,4,13,8,12,9,1,5,11,2,7,6,3,15),]
colnames(table_s5) <- c('','ClinGen HI genes','MHKO genes','SFARI ASD genes','DDG2P genes','eyeG2P genes','skinG2P genes')

table_s6_data <- auroc_test[[1]][c(1,3,5,9,6,7),]
table_s6_data1 <- auroc_test[[1]][c(1,3,5,9,6,7),][,c(1:14)]
table_s6_data1 <- data.frame(apply(table_s6_data1, 2, function(x) as.numeric(as.character(x))))
table_s6_data2 <- auroc_test[[3]][c(1,3,5,9,6,7),]
table_s6 <- c(colnames(table_s6_data2),'p-value')
for(i in 1:6){
  tmp <- c(paste(format(table_s6_data1[i,],nsmall=3),'±',format(table_s6_data2[i,],nsmall=3),sep = ''),table_s6_data[i,17])
  table_s6 <- rbind(table_s6,tmp)  
}
table_s6 <- as.data.frame(t(table_s6))
table_s6 <- table_s6[c(14,10,4,13,8,12,9,1,5,11,2,7,6,3,15),]
colnames(table_s6) <- c('','ClinGen HI genes','MHKO genes','SFARI ASD genes','DDG2P genes','eyeG2P genes','skinG2P genes')

write.table(table_s5,'table_s5.csv',quote = F,row.names = F,sep=',')
write.table(table_s6,'table_s6.csv',quote = F,row.names = F,sep=',')

######################
#Supplementary table 7
######################

annotation <- fread( "./other_data/anno_for_gene.txt")
gnomad <- fread('./other_data/gnomad.v2.1.1.lof_metrics.by_gene.txt')[,c('gene_id','exp_lof')]

# 988 tot  758 available
all_hi <- meta1[meta1$hi==1|meta1$asd==1|meta1$mhko==1|meta1$eyeg2p==1|meta1$sking2p==1|meta1$ddg2p==1,]
all_hi <- all_hi[!(is.na(all_hi$deeplof)&is.na(all_hi$oe_lof_upper)),]

all_hs <- meta1[meta1$hs==1,]

# low constrain but likely HI   352 tot  278 available
low_constrain <- all_hi[all_hi$oe_lof_upper>0.35 & all_hi$deeplof<0.835,]  

#short   135 tot  73 available
short <- all_hi[all_hi$gene_id %in% gnomad$gene_id[gnomad$exp_lof<10],]

#high constrain but unlikely HI  273 tot 209 available
high_constrain <- all_hs[all_hs$oe_lof_upper<0.35 & all_hs$deeplof>0.835,]  
high_constrain <- high_constrain[high_constrain$ratio2==0,]
hpo <- fread('./other_data/hpo_data.txt')
high_constrain <- high_constrain[!high_constrain$hgnc_symbol %in% hpo$V2,]

hl_count <- data.frame()
all_method_names <- colnames(all_hi)[17:30]
for (i in all_method_names){
  if(i=='oe_lof_upper'|i=='rvis'|i=='oe_mis_upper'|i=='virlof'){
    top <- meta1[order(meta1[i]),][1:2402,]
    top <- top$gene_id
    low_count <- sum(low_constrain$gene_id %in% top,na.rm = T)
    high_count <- sum(high_constrain$gene_id %in% top, na.rm = T )
    short_count <- sum(short$gene_id %in% top,na.rm = T)
    low_count <- c(i,low_count,'low constrain, likely HI')
    high_count <- c(i,high_count,'high constrain,unlikely HI')
    short_count <- c(i,short_count,'short, likely HI')
    count <- rbind.data.frame(low_count,high_count,short_count)
    colnames(count) <- c('methods','counts','conditions')
    hl_count <- rbind.data.frame(hl_count,count)
  }
  else{
    top <- meta1[order(-meta1[i]),][1:2402,]
    top <- top$gene_id
    low_count <- sum(low_constrain$gene_id %in% top,na.rm = T)
    high_count <- sum(high_constrain$gene_id %in% top, na.rm = T )
    short_count <- sum(short$gene_id %in% top,na.rm = T)
    low_count <- c(i,low_count,'low constrain, likely HI')
    high_count <- c(i,high_count,'high constrain,unlikely HI')
    short_count <- c(i,short_count,'short, likely HI')
    count <- rbind.data.frame(low_count,high_count,short_count)
    colnames(count) <- c('methods','counts','conditions')
    hl_count <- rbind.data.frame(hl_count,count)
  }
}

hl_count$conditions <- factor(hl_count$conditions,levels = c('short, likely HI','low constrain, likely HI','high constrain,unlikely HI'))
hl_count <- hl_count[order(hl_count$conditions),]
hl_count$counts <- as.numeric(hl_count$counts)
hl_count$methods <- rep(c('LOEUF','Mis. OEUF','CDS length','HIS','IS','GHIS','RVIS','HIPred','VIRLoF','Episcore','EDS','pHaplo','DeepLOF','DosaCNV-HI'),3)

short_data <- hl_count[hl_count$conditions=='short, likely HI',]
low_data <-  hl_count[hl_count$conditions=='low constrain, likely HI',]
high_data <-  hl_count[hl_count$conditions=='high constrain,unlikely HI',]

short_data <- cbind.data.frame(short_data, t(apply(short_data,1, function(p) prop.test(as.numeric(unlist(p[2])),nrow(short))$conf.int)))
short_data$inclusion_percent <- short_data$counts/nrow(short)
low_data <- cbind.data.frame(low_data, t(apply(low_data,1, function(p) prop.test(as.numeric(unlist(p[2])),nrow(low_constrain))$conf.int)))
low_data$inclusion_percent <- low_data$counts/nrow(low_constrain)
high_data <- cbind.data.frame(high_data, t(apply(high_data,1, function(p) prop.test(as.numeric(unlist(p[2])),nrow(high_constrain))$conf.int)))
high_data$inclusion_percent <- high_data$counts/nrow(high_constrain)

colnames(short_data)[4:5] <- c('lower','upper')
colnames(low_data)[4:5] <- c('lower','upper')
colnames(high_data)[4:5] <- c('lower','upper')

short_data$methods <- factor(short_data$methods,levels = short_data$methods[order(short_data$inclusion_percent)])
high_data$methods <- factor(high_data$methods,levels = high_data$methods[order(-high_data$inclusion_percent)])
low_data$methods <- factor(low_data$methods,levels = low_data$methods[order(low_data$inclusion_percent)])

short_data$se <- (short_data$upper-short_data$inclusion_percent)/1.96
high_data$se <- (high_data$upper-high_data$inclusion_percent)/1.96
low_data$se <- (low_data$upper-low_data$inclusion_percent)/1.96

short_data$short_inclusion_percent <- paste(format(round(short_data$inclusion_percent,3),nsmall = 3),'±',
                                            format(round(short_data$se,3),nsmall = 3), sep = '')

low_data$low_inclusion_percent <- paste(format(round(low_data$inclusion_percent,3),nsmall = 3),'±',
                                            format(round(low_data$se,3),nsmall = 3), sep = '')

high_data$high_inclusion_percent <- paste(format(round(high_data$inclusion_percent,3),nsmall = 3),'±',
                                            format(round(high_data$se,3),nsmall = 3), sep = '')

table_s7 <- cbind.data.frame(short_data$methods,short_data$short_inclusion_percent,low_data$low_inclusion_percent,high_data$high_inclusion_percent)
table_s7 <- table_s7[c(13,14,4,10,8,2,12,11,5,9,6,1,7,3),]
colnames(table_s7) <- c('methods','short','low','high')
short_p <- prop.test(c(29,34),c(73,73),correct = F)$p.value
low_p <- prop.test(c(122,115),c(278,278),correct = F)$p.value
high_p <- prop.test(c(31,15),c(209,209),correct = F)$p.value
p_row <- data.frame(methods = "p-value", short = short_p, low = low_p, high = high_p)
table_s7 <- rbind.data.frame(table_s7,p_row)
colnames(table_s7) <- c('','Short_likely HI genes','Low constraint_likely HI genes','High constraint_unlikely HI genes')

write.table(table_s7,'table_s7.csv',quote = F,row.names = F,sep=',')
