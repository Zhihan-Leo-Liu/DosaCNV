library(stringr)
library(data.table)
library(pROC)
library(readxl)
library(PRROC)
library(ggplot2)
library(coin)
library(grid)
library(cowplot)
library(egg)
library(MatchIt)

######################
# Figure 4a,b,c,d,e,f: AUC values for DosaCNV and 12 alternative methods across six HI and likely gene sets.
# Figure 4g,h,i: Inclusion percentage for DosaCNV and 12 alternative methods across three challenging gene sets.
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
    upper <- average + (upper-average)/1.96
    
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

# plot function for meta auroc
meta_auroc_plot <- function(data,target){
  a <- data
  method_names <- names(a[[2]])
  b <- as.data.frame(t(rbind.data.frame(a[[2]][target,], 
                                        a[[1]][target,][-c(length(a[[1]][target,]),length(a[[1]][target,])-1,length(a[[1]][target,])-2)],
                                        a[[3]][target,])))
  b[,1:3] <- apply(b[,1:3],2,as.numeric)
  colnames(b) <- c('lower','average','upper')
  b$methods <- method_names
  b
  b <- transform(b,methods=reorder(methods,average))
  str(b)
  
  meta_auroc_plot <- ggplot(b,aes(x=methods,y=average)) +
    geom_point(position=position_dodge(w=0.5),size=0.8,shape=15) + geom_errorbar(aes(ymin=lower,ymax=upper),size=0.3,width=0,position=position_dodge(.9))+ 
    ylim(0.42,1) + theme_classic() + ylab('AUC') + xlab(element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size=6),
          plot.title = element_text(face='bold',hjust=0.5),
          plot.margin=margin(0.2,0.2,0.2,0.2, unit="cm")) + coord_flip()
  return(meta_auroc_plot)
}

fig4a <- meta_auroc_plot(auroc_all,'hi') +  ggtitle(expression(paste(bold('ClinGen HI genes'))))  + theme(plot.title = element_text(size=6))
fig4b <-  meta_auroc_plot(auroc_all,'mhko') + ggtitle(expression(paste(bold('MHKO genes'))))  + theme(plot.title = element_text(size=6))
fig4c <-  meta_auroc_plot(auroc_all,'asd') + ggtitle(expression(paste(bold('SFARI ASD genes'))))  + theme(plot.title = element_text(size=6))
fig4d <-  meta_auroc_plot(auroc_all,'ddg2p') + ggtitle(expression(paste(bold('DDG2P genes')))) + theme(plot.title = element_text(size=6))
fig4e <- meta_auroc_plot(auroc_all,'eyeg2p') + ggtitle(expression(paste(bold('eyeG2P genes')))) + theme(plot.title = element_text(size=6))
fig4f <- meta_auroc_plot(auroc_all,'sking2p') + ggtitle(expression(paste(bold('skinG2P genes')))) + theme(plot.title = element_text(size=6))

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


fig4g <- ggplot(short_data,aes(x=methods,y=inclusion_percent)) + 
  geom_errorbar(aes(ymin=lower,ymax=inclusion_percent+se),width=0.3,size=0.35,position=position_dodge(w=0)) + 
  geom_bar(stat="identity", position=position_dodge(),fill="gray45",width = 0.8) + theme_classic() +
  coord_flip()  + scale_y_continuous(limits=c(0,1),expand = c(0,0))+
  #geom_text(aes(label=format(round(inclusion_percent, 3), nsmall = 3)), size=1.5,hjust=1.1, color="white", position = position_dodge(0.9)) +
  theme(text = element_text(size = 6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face='bold',hjust = 0.5),
        legend.title = element_blank(),
        plot.margin=margin(0.2,0.2,0.2,0.2, unit="cm")) + xlab(element_blank()) +
  ylab('Inclusion percentage')+
  ggtitle(expression(paste(bold("Short, likely HI genes"),' (n=73)')))  + theme(plot.title = element_text(size=6))

fig4h <- ggplot(low_data,aes(x=methods,y=inclusion_percent)) + 
  geom_errorbar(aes(ymin=lower,ymax=inclusion_percent+se),width=0.3,size=0.35,position=position_dodge(w=0)) + 
  geom_bar(stat="identity", position=position_dodge(),fill="gray45",width = 0.8) + theme_classic() +
  coord_flip()  + scale_y_continuous(limits=c(0,1),expand = c(0,0)) +
  #geom_text(aes(label=format(round(inclusion_percent, 3), nsmall = 3)), size=1.5,hjust=1.1, color="white", position = position_dodge(0.9)) +
  theme(text = element_text(size = 6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face='bold',hjust = 0.5),
        legend.title = element_blank(),
        plot.margin=margin(0.2,0.2,0.2,0.2, unit="cm")) + xlab(element_blank()) +
  ylab('Inclusion percentage')+
  ggtitle(expression(paste(bold("Low constraint, likely HI genes"),' (n=278)')))  + theme(plot.title = element_text(size=6))

fig4i <- ggplot(high_data,aes(x=methods,y=inclusion_percent)) + 
  geom_errorbar(aes(ymin=lower,ymax=inclusion_percent+se),width=0.3,size=0.35,position=position_dodge(w=0)) + 
  geom_bar(stat="identity", position=position_dodge(),fill="gray45",width = 0.8) + theme_classic() +
  coord_flip()  + scale_y_continuous(limits=c(0,1),expand = c(0,0)) +
  #geom_text(aes(label=format(round(inclusion_percent, 3),nsmall = 3)), size=1.5, hjust=1.1, color="white", position = position_dodge(0.9)) +
  theme(text = element_text(size = 6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face='bold',hjust = 0.5),
        legend.title = element_blank(),
        plot.margin=margin(0.2,0.2,0.2,0.2, unit="cm")) + xlab(element_blank()) +
  ylab('Inclusion percentage')+
  ggtitle(expression(paste(bold('High constraint, unikely HI genes'),' (n=209)')))  + theme(plot.title = element_text(size=6))

fig4 <- plot_grid(fig4a,fig4b,fig4c,
                  fig4d,fig4e,fig4f,
                  fig4g,fig4h,fig4i,
                  nrow = 3,
                  labels = 'AUTO',
                  label_size = 8,
                  rel_widths = c(1),
                  rel_heights = c(1))

pdf(file = './fig4.pdf',   # The directory you want to save the file in
    width = 7.2, # The width of the plot in inches
    height = 6) # The height of the plot in inches
fig4
dev.off()

#################
# Figure 5b: percentile plot + SHAP values heat map for 30 genes overlapped with 5 microdeletion
#################

our_score <- merge(our,meta[,c("gene_id",'hgnc_symbol')],by='gene_id')
our_score <- our_score[order(-our_score$DosaCNV_HI_score),]
our_score$rank <- 1:nrow(our_score)
our_score$percentile <- our_score$rank/nrow(our_score)

genes_1q21.1 <- c('BCL9','PRKAB2','GJA5','CHD1L','GJA8','ACP6','FMO5')
genes_15q13.3 <- c('FAN1','MTMR10','TRPM1','KLF13','OTUD7A','CHRNA7')
genes_17p12_hnpp_cmt1a <- c('COX10','CDRT15','HS3ST3B1','PMP22','TEKT3','CDRT4','TVP23C-CDRT4','TVP23C')
genes_17q21.31 <- c('CRHR1','SPPL2C','MAPT','STH','KANSL1')
genes_22q13<- c('MAPK8IP2','ARSA','SHANK3','ACR')

all_genes <- c(genes_1q21.1, genes_15q13.3, genes_17p12_hnpp_cmt1a, genes_17q21.31, genes_22q13)

# Extract subset of the data
selected_genes_data <- our_score[our_score$hgnc_symbol %in% all_genes,]
gene_percentile_df <- our_score[our_score$gene_id %in% selected_genes_data$gene_id, c('hgnc_symbol', 'percentile')]
gene_percentile_df$hgnc_symbol <- factor(gene_percentile_df$hgnc_symbol, levels = all_genes)

percentile_plot <- ggplot(gene_percentile_df, aes(x=hgnc_symbol, y=0, fill = percentile)) +
  geom_tile(colour='white', size=0.0)  + theme_bw() + xlab('') + ylab('') +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colors = c('#7a0177','#c51b8a','#f768a1','#f768a1','#fbb4b9','#fbb4b9','#feebe2','#feebe2'),
                       limits=c(0,1)) + theme(axis.text.x.top = element_text(vjust = 0.5)) +
  theme(legend.text=element_text(face="bold"),axis.text.x = element_text(size=11,angle = 90, vjust = 0.5, hjust=0.0),
        axis.ticks=element_line(size=0.4), axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        plot.background=element_blank(), aspect.ratio = 0.033) + 
  labs(fill = "DosaCNV-HI \n percentile")

y <- 7.8
y1 <- y+0.8

percentile_plot <- percentile_plot +  annotation_custom(grob = grid::textGrob(label = "1q21.1", hjust=0.3, gp=gpar(cex=1.2)),
                             xmin = 3.7, xmax = 3.7, ymin = y1, ymax = y1)  + coord_cartesian(clip="off") + 
  annotation_custom(grob = linesGrob(gp=gpar(lwd=2)), xmin = 0.65, xmax = 7.5, ymin = y, ymax = y) +
  
  annotation_custom(grob = grid::textGrob(label = "15q13.3", hjust=0.3, gp=gpar(cex=1.2)),
                    xmin = 10, xmax = 10, ymin = y1, ymax = y1) +
  annotation_custom(grob = linesGrob(gp=gpar(lwd=2)), xmin = 7.8, xmax = 13.5, ymin = y, ymax = y) +
  
  annotation_custom(grob = grid::textGrob(label = "17p12", hjust=0.3, gp=gpar(cex=1.2)),
                    xmin = 17.2, xmax = 17.2, ymin = y1, ymax = y1) +
  annotation_custom(grob = linesGrob(gp=gpar(lwd=2)), xmin = 13.8, xmax = 21.6, ymin = y, ymax = y) +
  
  annotation_custom(grob = grid::textGrob(label = "17q21.31", hjust=0.3, gp=gpar(cex=1.2)),
                    xmin = 23.4, xmax = 23.4, ymin = y1, ymax = y1) +
  annotation_custom(grob = linesGrob(gp=gpar(lwd=2)), xmin = 21.9, xmax = 26.5, ymin = y, ymax = y) +
  
  annotation_custom(grob = grid::textGrob(label = "22q13", hjust=0.3, gp=gpar(cex=1.2)),
                    xmin = 28.2, xmax = 28.2, ymin = y1, ymax = y1) +
  annotation_custom(grob = linesGrob(gp=gpar(lwd=2)), xmin = 26.8, xmax = 30.5, ymin = y, ymax = y)

shap <- fread('./other_data/all_shap.csv')
colnames(shap)[5:32] <- gsub('_',' ',colnames(shap)[5:32])
colnames(shap)[c(19,24,27)] <- c('GO embryo dev.','Tissue specificity (tau)','H2A.Z')

shap <- merge(shap,meta[,c("gene_id",'hgnc_symbol')],by='gene_id',all.x = T)
shap <- shap[!duplicated(shap),]

shap_df <- shap[shap$hgnc_symbol %in% gene_percentile_df$hgnc_symbol,-1]
shap_df_long <- melt(shap_df,id.vars = 'hgnc_symbol')
colnames(shap_df_long) <- c('hgnc_symbol','features','value')

feat_rank <- c('GTEx mean expression', 'Num tissues expressed',
               'CRISPR HAP1 KO qvalue', 'Num paralogs', 'DANN score',
               'PPI degree', 'Num interaction partners', 'mis_z', 'H3K4me3',
               'Promoter CpG density', 'Exonic phastCons100way', 'sHet',
               'Mis.OEUF', 'H3K27me3', 'Promoter phastCons100way', 'H2A.Z',
               'Num coding exons', 'Percent pathogenic SNVs', 'GO embryo dev.',
               'Num enhancers', 'Genic phastCons30way', 'lof_z', 'H3K9ac',
               'CADD raw', 'Genic phastCons100way', 'Tissue specificity (tau)',
               'LOEUF', 'MetaLR score', 'pLI', 'Num transcript isoforms',
               'Presence pathogenic SNVs')


shap_df_long$features <- factor(shap_df_long$features, levels = feat_rank)
shap_df_long$hgnc_symbol <- factor(shap_df_long$hgnc_symbol,levels = all_genes)
color <- c('#7a0177','#c51b8a','#f768a1','#fbb4b9','#feebe2',
           '#eff3ff','#bdd7e7','#6baed6','#3182bd','#08519c')

shap_plot <- ggplot(shap_df_long, aes(x=hgnc_symbol,y=features,fill=value)) + theme_bw()+
  geom_tile(colour='white',size=0.0) + coord_fixed(ratio=1) + xlab('') + ylab('') +
  scale_fill_gradientn(colors = rev(color),limits=c(-0.2,0.2),na.value ='#7a0177'  )+
  theme(axis.ticks.x =   element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_text(size=11),
        plot.background=element_blank()) + labs(fill = "SHAP value")


fig5b <- ggarrange(percentile_plot, shap_plot, 
              ncol = 1,
              heights = c(1,1))

pdf(file = './fig5b.pdf',   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
fig5b
dev.off()
