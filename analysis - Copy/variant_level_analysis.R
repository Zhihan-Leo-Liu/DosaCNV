library(data.table)
library(MatchIt)
library(pROC)
library(ggplot2)
library(reshape)
library(scales)
library(cowplot)
library(dplyr)

#############
# Figure 2A,C: ROC plots for unmatched and matched both (num of affected genes and deletion length) ClinVar held-out test
# set across methods.
# Figure 2B,D: Dot plots for AUC values of four length groups across methods.
#############

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

fig2a <- auc_plot(var_comp) + 
  ggtitle('Unmatched') + 
  theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))

fig2c <- auc_plot(var_comp_match_both) + 
  ggtitle('Matched') + 
  theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))

fig2b <- size_group_plot(size_group_comp(combined_variant,'unmatch')) + 
  ggtitle('Unmatched') + 
  theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))

fig2d <- size_group_plot(size_group_comp(combined_variant,'both')) + 
  ggtitle('Matched') + 
  theme(plot.title = element_text(size=7,face='bold',hjust = 0.5,vjust=-0.2))

fig2 <- plot_grid(fig2a,fig2b,
                  fig2c,fig2d,
                  nrow = 2,
                  labels = 'AUTO',
                  label_size = 8,
                  rel_widths = c(0.7,1))

fig2

pdf(file = './fig2.pdf',   # The directory you want to save the file in
    width = 7.2, # The width of the plot in inches
    height = 6) # The height of the plot in inches
fig2
dev.off()

####################################
# Figure 3A,B: Enrichment of case-associated deletions in each method's top predictions.
# Figure 3C: Scatter plot of DosaCNV scores for a set of Clinically curated case-only deletions. 
####################################

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

nstd100_plot <- ggplot(nstd100_data, aes(x = method, y = exp(or))) + 
  geom_errorbar(aes(ymin = exp(or) - 0.5, ymax = upper), width = 0.3, size = 0.3, position = position_dodge(0)) + 
  geom_bar(stat = "identity", position = position_dodge(), fill = "gray45", width = 0.8) +
  coord_flip() + 
  scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 6, face = 'bold', hjust = 1), legend.title = element_blank(), plot.margin = margin(0.3, 0.3, 0.3, 0, unit = "cm")) +
  xlab(element_blank()) + ylab('Enrichment (odds ratio)') + ggtitle(expression(paste(bold("nstd100: Case-control deletions"), " (Coe et al., ", italic("Nat Genet"), ", 2014)")))

nstd173_plot <- ggplot(nstd173_data, aes(x = method, y = exp(or) + 0.005)) + 
  geom_errorbar(aes(ymin = exp(or), ymax = upper), width = 0.3, size = 0.3, position = position_dodge(0)) + 
  geom_bar(stat = "identity", position = position_dodge(), fill = "gray45", width = 0.8) +
  coord_flip() + 
  scale_y_continuous(trans = log2_trans(), breaks = c(1,2^0.33,2^0.66,2^1), labels = trans_format("log2", math_format(2^.x))) +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 6, face = 'bold', hjust = 1), legend.title = element_blank(), plot.margin = margin(0.3, 0.3, 0.3, 0, unit = "cm")) +
  xlab(element_blank()) + ylab('Enrichment (odds ratio)') + ggtitle(expression(paste(bold("nstd173: Case-control deletions"), " (Zarrei et al., ", italic("NPJ Genom Med"), ", 2019)")))

fig3a <- plot_grid(nstd100_plot,nstd173_plot,
                   nrow = 2,
                   labels = 'AUTO',
                   label_size = 8)

fig3a

######################################################## 
########################################################

# Read and filter data
sickkids_data <- fread('./other_data/sickkid_data.bed')
sickkids_data <- sickkids_data[sickkids_data$V4 == 'DEL',]
sickkids_data$V5 <- paste('new_', seq(1:nrow(sickkids_data)), sep = '')
colnames(sickkids_data) <- c('chr', 'start', 'end', 'type', 'variant_id', 'label')

# Merge with variant scores
sickkid_score <- fread('./our_scores/sickkid_score.tsv')
colnames(sickkid_score)[2] <- 'our_score'
merged_data <- merge(sickkid_score, sickkids_data, by = 'variant_id')

# Convert chromosome names to numbers
merged_data$chr_num <- as.numeric(sub('chr', '', merged_data$chr))

# Create chromosome length dataframe
chr_lengths <- data.frame(
  chr_length = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432),
  chr = paste('chr', 1:22, sep = '')
)

# Merge data with chromosome lengths
merged_data <- merge(merged_data, chr_lengths, by = 'chr')

# Calculate normalized chromosome positions
# Calculate minimum start positions for each chromosome
min_positions <- aggregate(start ~ chr_num, data = merged_data, FUN = min)
colnames(min_positions) <- c("chr_num", "min")

# Calculate maximum start positions for each chromosome
max_positions <- aggregate(start ~ chr_num, data = merged_data, FUN = max)
colnames(max_positions) <- c("chr_num", "max")

# Merge the min and max values together
norm_positions <- merge(min_positions, max_positions, by = "chr_num")

# Merge these normalized positions back with the main dataset
merged_data <- merge(merged_data, norm_positions, by = "chr_num", all.x = TRUE)

merged_data$location <- ((merged_data$start - merged_data$min) / (merged_data$max + 1 - merged_data$min)) + merged_data$chr_num

# Adjust labels for clarity
merged_data$label[merged_data$label == 'B'] <- 'Benign'
merged_data$label[merged_data$label == 'H'] <- 'Pathogenic'

# Classify data based on our_score
classified_data <- merged_data
classified_data$label_hat[classified_data$our_score >= 0.8] <- 1
classified_data$label_hat[classified_data$our_score < 0.8] <- 0

# Create confusion matrix  TPR=0.970(32/33)  FPR = 0.028(21/758) Precision = 0.604(32/53)
confusion_matrix <- table(classified_data$label, classified_data$label_hat)
roc_object <- roc(merged_data$label,merged_data$our_score)  # AUC = 0.9981

fig3b <- ggplot(merged_data, aes(x = location, y = our_score, col = label)) +
  geom_point(position = position_jitter(h = 0.005, w = 0.05), size = 1, shape = 1) + 
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = 'dashed') +
  scale_color_manual(values = c('gray34', "red")) +
  scale_x_discrete(limits = 1:22) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = expansion(add = 0.05)) +
  ggtitle(expression(paste(bold("Hospital of Sick Children: Case deletions"), " (Foong et al., ", italic("PLoS One"), ", 2015)"))) +
  theme(
    text = element_text(size = 6),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 6, face = 'bold', hjust = 0.5),
    axis.text.x = element_text(size = 6, hjust = -0.5),
    legend.position = c(.93, 0.6),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "cm"),
    legend.title = element_blank(),
    legend.key.size = unit(0.15, 'cm'),
    legend.background = element_rect(fill = "white", color = "black")
  ) + 
  ylab('DosaCNV score') + 
  xlab('Chromosome')


fig3b

fig3 <- plot_grid(fig3a,fig3b,
                  nrow = 1,
                  labels = c('','C'),
                  label_size = 8,
                  rel_widths = c(0.9,1))
fig3

pdf(file = './fig3.pdf',   # The directory you want to save the file in
    width = 7.2, # The width of the plot in inches
    height = 4) # The height of the plot in inches
fig3
dev.off()


a <- roc(merged_data$label,merged_data$our_score)
which(a$thresholds==0.8)

###################################
# Enrichment of true pathogenic deletions in DosaCNV's top predictions for two case-control datasets(nstd100, nstd173).
###################################

nstd100_102_overlap <- as.data.frame(fread('./overlap_analysis/nstd100_102_overlap.bed'))[,c(4,8,9)]
nstd173_102_overlap <- as.data.frame(fread('./overlap_analysis/nstd173_102_overlap.bed'))[,c(4,8,9)]

colnames(nstd100_102_overlap) <- c('nstd100_id','nstd102_id','label')
colnames(nstd173_102_overlap) <- c('nstd173_id','nstd102_id','label')

nstd100_102_uniq <- nstd100_102_overlap %>%
  group_by(nstd100_id) %>%
  summarise(avg = mean(label))

nstd173_102_uniq <- nstd173_102_overlap %>%
  group_by(nstd173_id) %>%
  summarise(avg = mean(label))

top5_nstd100 <- our_nstd100_score[order(-our_nstd100_score$our_score),][1:1454,]
top5_nstd173 <- our_nstd173_score[order(-our_nstd173_score$our_score),][1:134,]

nstd100_102_patho <- nstd100_102_uniq[nstd100_102_uniq$avg==1,]           # 164
nstd100_102_benign <- nstd100_102_uniq[nstd100_102_uniq$avg==0,]         # 172 
nstd173_102_patho <- nstd173_102_uniq[nstd173_102_uniq$avg==1,]           # 13
nstd173_102_benign <- nstd173_102_uniq[nstd173_102_uniq$avg==0,]         # 52

summary(nstd100_102_patho$nstd100_id %in% top5_nstd100$variant_id)   # 148/164 in top5 nstd100 predictions of dosacnv  
summary(nstd173_102_patho$nstd173_id %in% top5_nstd173$variant_id)   # 4/13 in top5 nstd173 predictions of dosacnv
summary(nstd100_102_benign$nstd100_id %in% top5_nstd100$variant_id)  # 0/168  '' ''
summary(nstd173_102_benign$nstd173_id %in% top5_nstd173$variant_id)  # 0/49   '' ''

nstd100_case_id <- meta_nstd100_score[meta_nstd100_score$status==1,'variant_id']
nstd173_case_id <- meta_nstd173_score[meta_nstd173_score$status==1,'variant_id']

summary(nstd100_102_patho$nstd100_id %in% nstd100_case_id$variant_id)   # 154/164 in cases
summary(nstd173_102_patho$nstd173_id %in% nstd173_case_id$variant_id)   # 11/13 in cases
