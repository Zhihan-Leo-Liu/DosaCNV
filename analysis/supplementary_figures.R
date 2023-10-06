library(data.table)
library(MatchIt)
library(pROC)
library(ggplot2)
library(reshape)
library(scales)
library(cowplot)

#######################
#Supplementary figure 1
#######################

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

our_nstd102_score <- fread('./our_scores/nstd102_score.tsv')
meta_nstd102_scores <- fread('./meta_scores/meta_nstd102_scores.csv')

colnames(our_nstd102_score)[2] <- 'our_score'
meta_nstd102_scores$var_length <- meta_nstd102_scores$end_hg19 - meta_nstd102_scores$start_hg19

combined_variant <- merge(our_nstd102_score, meta_nstd102_scores, by = 'variant_id', all.x = T, all.y = F)
combined_variant <- combined_variant[complete.cases(combined_variant), ]
combined_variant <- combined_variant[, c(2,11,12,13,14,15,17,19,18)]

# Match Data by Different Criteria
set.seed(1)

# Match by gene
matched_gene <- match.data(matchit(label ~ num_coding_gene, data = combined_variant, method = 'nearest', caliper = 0, discard = 'both', m.order = 'random'))
matched_gene <- matched_gene[, c(7,8,9)]
matched_gene$matching <- 'Matched by num of affected genes'

# Match by deletion length
matched_length <- match.data(matchit(label ~ var_length, data = combined_variant, method = 'nearest', caliper = 1, discard = 'both', m.order = 'random'))
matched_length <- matched_length[, c(7,8,9)]
matched_length$matching <- 'Matched by deletion length'

# Match by both gene and length
matched_both <- match.data(matchit(label ~ var_length, data = combined_variant, exact = ~ num_coding_gene, method = 'nearest', caliper = 1, discard = 'both', m.order = 'random'))
matched_both <- matched_both[, c(7,8,9)]
matched_both$matching <- 'Matched by both'

# Prepare Unmatched Data
unmatched <- combined_variant[, c(7,8,9)]
unmatched$matching <- 'Unmatched'

# Labeling
data_list <- list(unmatched, matched_gene, matched_length, matched_both)
data_list <- lapply(data_list, function(df) {
  df$label[df$label == 1] <- 'Pathogenic'
  df$label[df$label == 0] <- 'Benign'
  return(df)
})

# Combine All Data
dist_data <- do.call(rbind.data.frame, data_list)
dist_data$matching <- factor(dist_data$matching, levels = c('Unmatched', 'Matched by deletion length', 'Matched by num of affected genes', 'Matched by both'))

common_theme <- theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=7),
        legend.title = element_blank(),legend.key.size = unit(0.5,'cm'),
        legend.margin = margin(-8,-0.5,-0.5,-0.5), legend.spacing.x = unit(0.5,'mm'),legend.spacing.y = unit(0.5,'mm'),
        plot.margin = margin(0.3,0.3,0.3,0.3, unit="cm"))

figs1a <- ggplot(dist_data, aes(x=matching, y=log(num_coding_gene), fill=label)) + 
  geom_split_violin(trim=T, alpha=0.5) + 
  scale_fill_manual(values=c('Pathogenic'= "red", 'Benign'="gray34")) + 
  xlab('') + ylab('log num of affected genes') + 
  common_theme

figs1b <- ggplot(dist_data, aes(x=matching, y=log(var_length), fill=label)) + 
  geom_split_violin(trim=T, alpha=0.5) + 
  scale_fill_manual(values=c('Pathogenic'= "red", 'Benign'="gray34")) + 
  xlab('') + ylab('log deletion length') + 
  common_theme

figs1 <- plot_grid(figs1a, NULL, figs1b, nrow = 3, labels = c('A', '', 'B'), label_size = 8, rel_heights = c(1,0.1,1))

pdf(file = './figs1.pdf',   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches
figs1
dev.off()

#######################
#Supplementary figure 2
#######################

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
  
  # Rank methods based on AUC
  methods_rank <- names(auc_value[order(-auc_value)])
  
  # Compute p-value for our_score vs. tada
  pvalue <- roc.test(roc_obj[['our_score']], roc_obj[['xcnv']])$p.value
  
  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, round(pvalue, 4)), nrow = 1))
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results unless quiet is TRUE
  if (!quiet) {
    print(list_result)
  }
  
  return(list(roc_obj, list_result))
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
  
  # Rank methods based on AUC (this isn't used later, but kept for consistency)
  methods_rank <- names(auc_value[order(-auc_value)])
  
  # Compute p-value for our_score vs. tada
  pvalue <- roc.test(roc_obj[['our_score']], roc_obj[['tada']])$p.value
  
  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, round(pvalue, 4)), nrow = 1))
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results
  print(list_result)
  
  return(list(roc_obj, list_result))
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
  
  # Rank methods based on AUC (this isn't used later, but kept for consistency)
  methods_rank <- names(auc_value[order(-auc_value)])
  
  # Compute p-value for our_score vs. tada
  pvalue <- roc.test(roc_obj[['our_score']], roc_obj[['tada']])$p.value
  
  # Combine results into a dataframe
  list_result <- as.data.frame(matrix(c(auc_value, round(pvalue, 4)), nrow = 1))
  colnames(list_result) <- c(methods_name, 'pvalue')
  
  # Display the results
  print(list_result)
  
  return(list(roc_obj, list_result))
}

auc_plot <- function(comparison) {
  
  # Extract ROC object from the comparison function
  roc_objects <- comparison(combined_variant)[[1]]
  
  # Define method names and initialize an empty list for AUCs
  method_names <- c('DosaCNV', 'CADD-SV', 'TADA', 'X-CNV', 'StrVCTVRE', 'HIS-LOD', 'num affected genes', 'deletion length')
  auc_list <- c()
  
  # Populate the AUC list for each method
  for (i in 1:8) {
    auc <- format(round(roc_objects[[i]]$auc, 3), nsmall = 3)
    auc_list <- c(auc_list, paste('(AUC = ', auc, ')', sep = ''))
  }
  
  # Combine method names with their AUCs
  methods_auc <- paste(method_names, auc_list)
  
  # Rename and reorder roc_objects based on methods_auc
  roc_objects <- setNames(roc_objects, methods_auc)
  roc_objects <- roc_objects[c(1, 3, 5, 6, 4, 2, 7, 8)]
  
  # Construct the ROC plot
  plot <- ggroc(roc_objects, legacy.axes = TRUE, size = 0.5) + 
    theme_classic() +
    scale_x_continuous(n.breaks = 6) +
    scale_y_continuous(n.breaks = 6) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75, 0.23),
      text = element_text(size = 7),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, 'cm'),
      legend.margin = margin(-8, -0.5, -0.5, -0.5),
      legend.spacing.x = unit(0, 'mm'),
      legend.spacing.y = unit(0.5, 'mm'),
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm")
    ) +
    scale_color_manual(values = c('black', 'coral', 'blueviolet', 'darkolivegreen3', 'deeppink2', 'royalblue3', 'deepskyblue2', 'darkgray')) +
    xlab('False positive rate') +
    ylab('True positive rate')
  
  return(plot)
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
  
  out_plot <- rbind.data.frame(q1_auc_full, q2_auc_full, q3_auc_full, q4_auc_full)
  out_plot$q <- rep(c('<50 ', '50-100 ', '100-1000 ', '>1000 '), each = nrow(out_plot) / 4)
  out_plot$method <- rep(colnames(scores)[1:8], 4)
  colnames(out_plot)[1:3] <- c('lower', 'average', 'upper')
  out_plot$upper <- out_plot$average + (out_plot$upper - out_plot$average) / 2 
  out_plot$lower <- out_plot$average - (out_plot$average - out_plot$lower) / 2
  
  num_q <- c(nrow(q1), nrow(q2), nrow(q3), nrow(q4))
  
  print(out)
  
  return(list(out, out_plot, num_q))
}

size_group_plot <- function(out) {
  
  outplot <- out[[2]]
  num_q <- out[[3]]
  
  # Method names and reformatting the 'q' variable
  method_names <- c('DosaCNV', 'CADD-SV', 'TADA', 'X-CNV', 'StrVCTVRE', 'HIS-LOD', 'num affected genes', 'deletion length')
  outplot$method <- rep(method_names, 4)
  
  outplot$q <- gsub(" ", "", paste(outplot$q, '(n='))
  outplot$q <- paste(outplot$q, rep(num_q, each = 8), ')', sep = '')
  
  # Ordering factors for plotting
  outplot$method <- factor(outplot$method, levels = c('DosaCNV', 'TADA', 'StrVCTVRE', 'HIS-LOD','X-CNV', 'CADD-SV','num affected genes', 'deletion length'))
  outplot$q <- factor(outplot$q, levels = unique(outplot$q))
  
  # Constructing the ggplot
  plot <- ggplot(outplot, aes(x = q, y = average, col = method)) + 
    geom_point(position = position_dodge(w = 0.5), size = 0.7, shape = 15) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.3, position = position_dodge(w = 0.5)) + 
    ylim(0.25, 1) + 
    theme_classic() + 
    ylab('AUROC') + 
    xlab(element_blank()) +
    scale_color_manual(values = c('DosaCNV' = 'black',
                                  'TADA' = 'coral',
                                  'StrVCTVRE' = 'blueviolet',
                                  'HIS-LOD' = 'darkolivegreen3',
                                  'X-CNV' = 'deeppink2',
                                  'CADD-SV' = 'royalblue3',
                                  'num affected genes' = 'deepskyblue2',
                                  'deletion length' = 'darkgrey')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size = 7), axis.text.x = element_text(size = 5.5),
          legend.title = element_blank(), legend.key.size = unit(0.3, 'cm'),
          legend.margin = margin(-8, -0.5, -0.5, -0.5), 
          legend.spacing.x = unit(0, 'mm'),
          legend.spacing.y = unit(0.5, 'mm'),
          plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm")) +
    xlab('Length (kb)') + 
    ylab('AUC') + 
    geom_hline(yintercept = 0.5, linetype = 2, size = 0.5)
  
  return(plot)
}

seed = 1
figs2a <- auc_plot(var_comp_match_cdg) + ggtitle('Matched by number of affected genes') + theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))
figs2c <- auc_plot(var_comp_match_length) + ggtitle('Matched by deletion length') + theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))
figs2b <- size_group_plot(size_group_comp(combined_variant,'cdg')) + ggtitle('Matched by number of affected genes') + theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))
figs2d <- size_group_plot(size_group_comp(combined_variant,'length')) + ggtitle('Matched by deletion length') + theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))
  
  
figs2 <- plot_grid(figs2a,figs2b,
                    figs2c,figs2d,
                    nrow = 2,
                    labels = 'AUTO',
                    label_size = 8,
                    rel_widths = c(0.7,1))

pdf(file = './figs2.pdf',   # The directory you want to save the file in
    width = 7.2, # The width of the plot in inches
    height = 6) # The height of the plot in inches
figs2
dev.off()

#######################
#Supplementary figure 3
#######################

meta_sickkid_scores <- as.data.frame(fread('./meta_scores/meta_sickkid_scores.csv'))
meta_sickkid_scores <- meta_sickkid_scores[,-c(1,2,3,4)]
meta_sickkid_scores$label <- as.numeric(meta_sickkid_scores$label)

methods_name <- colnames(meta_sickkid_scores)
colnames(meta_sickkid_scores) <- c('label','caddsv','tada','xcnv','strv','huang_score','num_coding_gene','del_length','our')

roc_object <- roc(label ~ ., data = meta_sickkid_scores)
auc_value <- lapply(roc_object, function(x) as.vector(ci(x)))

lower <- unlist(lapply(auc_value, function(l) l[[1]]))
auroc <- unlist(lapply(auc_value, function(l) l[[2]]))
upper <- unlist(lapply(auc_value, function(l) l[[3]]))

# Create result data frame
result <- data.frame(AUROC = auroc, lower = lower, upper = upper)
result$method <- methods_name[-1]
result$lower <- result$AUROC - (result$AUROC - result$lower) / 2
result$upper <- result$AUROC + (result$upper - result$AUROC) / 2
result$method <- factor(result$method, levels = result$method[order(result$AUROC)])

figs3 <- ggplot(result, aes(x = method, y = AUROC)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, size = 0.3, position = position_dodge(w = 0)) + 
  geom_bar(stat = "identity", position = position_dodge(), fill = "gray45", width = 0.8) + 
  coord_flip(ylim = c(0.8, 1)) + 
  geom_text(aes(label = format(round(AUROC, 3), nsmall = 3)), size = 1.5, hjust = 1.15, color = "white", position = position_dodge(0.9)) +
  theme_classic() +
  theme(
    text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 7, face = 'bold', hjust = 1.1),
    legend.title = element_blank(),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm")
  ) + 
  xlab(element_blank()) + 
  ylab('AUC') +
  ggtitle(expression(paste(bold("Hospital of Sick Children: Case deletions"), " (Foong et al., ", italic("PLoS One"), ", 2015)")))

roc.test(roc_object$our,roc_object$strv)

pdf(file = 'C:\\Users\\zhihan leo liu\\Desktop\\Huang\\Model\\multitask2\\figures\\figs3.pdf',   # The directory you want to save the file in
    width = 3.6, # The width of the plot in inches
    height = 2) # The height of the plot in inches
figs3
dev.off()

#######################
#Supplementary figure 4
#######################

meta_gene_scores <- fread('./meta_scores/meta_gene_scores.csv') 
meta_gene_scores[meta_gene_scores=='missing'] <- NA 
meta_gene_scores <- unique(as.data.frame(meta_gene_scores))
meta_gene_scores[, 3:ncol(meta_gene_scores)] <- lapply(meta_gene_scores[, 3:ncol(meta_gene_scores)], as.numeric)

our <- fread('./our_scores/precomputed_gene_score.tsv')[, c(1,2)]
meta_gene_scores <- merge(na.omit(meta_gene_scores), our, by ='gene_id', all.x=T)

colMeans(meta_gene_scores[meta_gene_scores$train==0,][,5:14])*nrow(meta_gene_scores[meta_gene_scores$train==0,])

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

auroc_test <- meta_auroc(meta_gene_scores[meta_gene_scores$train==0,])

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

figs4a <- meta_auroc_plot(auroc_test,'hi') +  ggtitle(expression(paste(bold('ClinGen HI genes'))))  + theme(plot.title = element_text(size=6))
figs4b <-  meta_auroc_plot(auroc_test,'mhko') + ggtitle(expression(paste(bold('MHKO genes'))))  + theme(plot.title = element_text(size=6))
figs4c <-  meta_auroc_plot(auroc_test,'asd') + ggtitle(expression(paste(bold('SFARI ASD genes'))))  + theme(plot.title = element_text(size=6))
figs4d <-  meta_auroc_plot(auroc_test,'ddg2p') + ggtitle(expression(paste(bold('DDG2P genes')))) + theme(plot.title = element_text(size=6))
figs4e <- meta_auroc_plot(auroc_test,'eyeg2p') + ggtitle(expression(paste(bold('eyeG2P genes')))) + theme(plot.title = element_text(size=6))
figs4f <- meta_auroc_plot(auroc_test,'sking2p') + ggtitle(expression(paste(bold('skinG2P genes')))) + theme(plot.title = element_text(size=6))


figs4 <- plot_grid(figs4a,figs4b,figs4c,
                  figs4d,figs4e,figs4f,
                  nrow = 2,
                  labels = 'AUTO',
                  label_size = 8,
                  rel_widths = c(1),
                  rel_heights = c(1))

pdf(file = './figs4.pdf',   # The directory you want to save the file in
    width = 7.2, # The width of the plot in inches
    height = 4) # The height of the plot in inches
figs4
dev.off()
