# load libraries
suppressPackageStartupMessages({
    library(caret)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(tidyverse)
    library(ggh4x)
    library(RColorBrewer)
})

set.seed(200)  


############################################################
# Load in data
############################################################

# load in stability data
load("../results/data/temp/gene_isoform_stability.RData")
load("../results/data/temp/circ_stability2.RData")

# load in circRNA features
load("../results/data/circ_stability_gcsi_features2.RData")
load("../results/data/circ_stability_ccle_features2.RData")
load("../results/data/circ_stability_gdsc_features2.RData")


############################################################
# Format circRNA dataframes
############################################################

# function to format circRNA stability dataframes
c_stability <- function(stability, gcsi, ccle, gdsc) {

  # reformat stability dataframe
  stability <- reshape2::dcast(stability, transcript ~ pair, value.var = "si")
  colnames(stability) <- c("transcript", "gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")
  
  # get median expression
  stability$gcsi_median <- stability$ccle_median <- stability$gdsc_median <- 0
  for (i in seq_along(stability$transcript)) {
    transcript <- stability$transcript[i]

    if (transcript %in% rownames(gcsi)) {
      stability$gcsi_median[i] <- gcsi$MedianExp[rownames(gcsi) == transcript]
    }
    if (transcript %in% rownames(ccle)) {
      stability$ccle_median[i] <- ccle$MedianExp[rownames(ccle) == transcript]
    }
    if (transcript %in% rownames(gdsc)) {
      stability$gdsc_median[i] <- gdsc$MedianExp[rownames(gdsc) == transcript]
    }
  }

  # get GC, nexons, and length
  feats <- data.frame(transcript = c(rownames(gcsi), rownames(ccle), rownames(gdsc)),
                      gc = c(gcsi$GC, ccle$GC, gdsc$GC),
                      nexons = c(gcsi$NExons, ccle$NExons, gdsc$NExons),
                      length = c(gcsi$Length, ccle$Length, gdsc$Length))
  feats <- unique(feats)
  feats <- feats[match(stability$transcript, feats$transcript),]
  
  stability$gc <- feats$gc
  stability$n_exon <- feats$nexons
  stability$length <- feats$length

  return(stability)
} 

ciri_stability <- c_stability(ciri_stability, ciri_gcsi_ft, ciri_ccle_ft, ciri_gdsc_ft)
circ_stability <- c_stability(circ_stability, circ_gcsi_ft, circ_ccle_ft, circ_gdsc_ft)
cfnd_stability <- c_stability(cfnd_stability, cfnd_gcsi_ft, cfnd_ccle_ft, cfnd_gdsc_ft)
fcrc_stability <- c_stability(fcrc_stability, fcrc_gcsi_ft, fcrc_ccle_ft, fcrc_gdsc_ft)

# save files
write.csv(ciri_stability, file = "../results/data/temp/ciri_stability.csv", quote = F, row.names = F)
write.csv(circ_stability, file = "../results/data/temp/circ_stability.csv", quote = F, row.names = F)
write.csv(cfnd_stability, file = "../results/data/temp/cfnd_stability.csv", quote = F, row.names = F)
write.csv(fcrc_stability, file = "../results/data/temp/fcrc_stability.csv", quote = F, row.names = F)


############################################################
# FROM HERE: RUN FeatureInfluence.ipynb
############################################################


############################################################
# Define function to load in multivariable model results
############################################################

# helper function to format each dataset pair result
load_pair <- function(path, pair, filetype) {

  # process model datafile
  if (filetype == "model") {
    
    # load in data files
    lm <- read.csv(paste0(path, pair, "lm.csv"))
    # remove brackets from Pearson
    lm$Pearson <- gsub("\\[", "", gsub("]", "", lm$Pearson)) |> as.numeric()
    lm$alpha <- lm$max_iter <- NA
    ls <- read.csv(paste0(path, pair, "lasso.csv"))

    # merge data files
    df <- rbind(lm, ls)
    df$pair <- gsub("/", "", pair)
  }

  # process feature datafile
  if (filetype == "features") {

    # load in data files
    lm <- read.csv(paste0(path, pair, "lm_features.csv"))
    ls <- read.csv(paste0(path, pair, "lasso_features.csv"))

    # add in labels
    lm$model <- "Linear"
    ls$model <- "LASSO"

    # merge data files
    df <- rbind(lm, ls)
    df$pair <- gsub("/", "", pair)

    # rename median
    df$Peak[df$Peak %in% c("gdsc_median", "ccle_median", "gcsi_median")] <- "median"
    colnames(df)[1] <- "Feature"
  }
  return(df)
}

# function to load in results
load_model <- function(label, filetype) {

  path = paste0("results/data/feature_influence/multivariable/", label)

  # process model files
  if (filetype == "model") {
    df <- rbind(
      load_pair(path, "/gcsi_ccle/", "model"),
      load_pair(path, "/gcsi_gdsc/", "model"),
      load_pair(path, "/gdsc_ccle/", "model")
    )
  }
  # process feature files
  if (filetype == "features") {
    df <- rbind(
      load_pair(path, "/gcsi_ccle/", "features"),
      load_pair(path, "/gcsi_gdsc/", "features"),
      load_pair(path, "/gdsc_ccle/", "features")
    )
  }
  df$label <- label
  return(df)
}


############################################################
# Load in mulivariable results
############################################################

# set up palette for plotting
pal = c("#DACCAB", "#C78B76", "#9D3737")

# set up model results
model_df <- rbind(
  load_model("Gene_Expression", "model"),
  load_model("Isoform_Expression", "model"),
  load_model("CIRI2", "model"),
  load_model("CIRCexplorer2", "model"),
  load_model("circRNA_finder", "model"),
  load_model("find_circ", "model")
)

# format model dataframe
model_df$Fold <- factor(model_df$Fold, levels = c(1, 2, 3, 4, 5))
model_df$pair <- factor(model_df$pair, 
  levels = c("gcsi_ccle", "gcsi_gdsc", "gdsc_ccle"),
  labels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))
model_df$label <- factor(model_df$label, 
  levels = c("Gene_Expression", "Isoform_Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"),
  labels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))


# set up feature results
feature_df <- rbind(
  load_model("Gene_Expression", "features"),
  load_model("Isoform_Expression", "features"),
  load_model("CIRI2", "features"),
  load_model("CIRCexplorer2", "features"),
  load_model("circRNA_finder", "features"),
  load_model("find_circ", "features")
)

# format feature dataframe
feature_df$Feature <- factor(feature_df$Feature, 
  levels = c("median", "gc", "n_exon", "length"),
  labels = c("Median Exp", "GC%", "No. Exons", "Length"))
feature_df$label <- factor(feature_df$label, 
  levels = c("Gene_Expression", "Isoform_Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"),
  labels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))


############################################################
# Format and save model results
############################################################

# function to compile results
compile_res <- function(model) {
  
  df <- model_df[model_df$Model == model, -c(5:6)]
  df$cat <- paste(df$label, df$pair, sep = "-")
  
  # get average and max
  res <- df %>%
    group_by(cat) %>%
    summarise(
      avg_spearman = mean(Spearman, na.rm = TRUE),
      sd_spearman  = sd(Spearman, na.rm = TRUE),
      avg_pearson  = mean(Pearson, na.rm = TRUE),
      sd_pearson   = sd(Pearson, na.rm = TRUE),
      max_spearman = Spearman[which.max(abs(Spearman))],
      max_pearson  = Pearson[which.max(abs(Pearson))],
      .groups = "drop"
    )
  
  # format results
  res$label <- gsub("-.*", "", res$cat)
  res$pair <- gsub(".*-", "", res$cat)
  res$label <- factor(res$label, levels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
  res$pair <- factor(res$pair, levels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))

  res <- res[order(res$label, res$pair),] |> as.data.frame()

  # transpose
  rownames(res) <- res$cat
  res <- t(res[,-c(1, 8, 9)]) |> as.data.frame()

  return(res)
}

# compile results
lm_compile <- compile_res("Linear")
ls_compile <- compile_res("LASSO")

# write results
write.csv(lm_compile, file = "results/data/feature_influence/multivariable/lm.csv", quote = F, row.names = T)
write.csv(ls_compile, file = "results/data/feature_influence/multivariable/LASSO.csv", quote = F, row.names = T)


############################################################
# Plot model results (across all folds)
############################################################

# function to plot model correlations
plot_model <- function(model, corr, scales = "fixed") {

  toPlot <- model_df[model_df$Model == model,]

  p <- ggplot(toPlot, aes(x = pair, y = .data[[corr]], fill = pair, alpha = Fold)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    scale_fill_manual(values = pal) +
    facet_wrap(.~label, nrow = 1, scales = scales) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm'),
          axis.text.x = element_text(angle = 90, hjust = 0.1)) +
    labs(fill = "Dataset Pair", x = "Dataset Pair") 

  png(paste0("results/figures/", model, "_model_", corr, "_", scales, "_scales.png"), width=11, height=5, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot model correlations
plot_model("Linear", "Spearman")
plot_model("Linear", "Spearman", "free_y")
plot_model("LASSO", "Spearman")
plot_model("LASSO", "Spearman", "free_y")

plot_model("Linear", "Pearson")
plot_model("Linear", "Pearson", "free_y")
plot_model("LASSO", "Pearson")
plot_model("LASSO", "Pearson", "free_y")


############################################################
# Plot model results (average across folds)
############################################################

# function to plot model correlations
plot_model <- function(df, model, corr, scales = "fixed") {

  toPlot <- t(df) |> as.data.frame()

  # format columns
  if (corr == "Spearman") {
    colnames(toPlot)[colnames(toPlot) == "avg_spearman"] <- "corr"
    colnames(toPlot)[colnames(toPlot) == "sd_spearman"] <- "sd"
  }
  if (corr == "Pearson") {
    colnames(toPlot)[colnames(toPlot) == "avg_pearson"] <- "corr"
    colnames(toPlot)[colnames(toPlot) == "sd_pearson"] <- "sd"
  }
  toPlot$corr <- as.numeric(toPlot$corr)
  toPlot$sd <- as.numeric(toPlot$sd)

  toPlot$pair <- factor(toPlot$pair, 
    levels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))
  toPlot$label <- factor(toPlot$label, 
    levels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))



  p <- ggplot(toPlot, aes(x = pair, y = corr, fill = pair)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=corr-sd, ymax=corr+sd), width=.2, position=position_dodge(.9)) +
    scale_fill_manual(values = pal) +
    facet_wrap(.~label, nrow = 1, scales = scales) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm'),
          axis.text.x = element_text(angle = 90, hjust = 0.1)) +
    labs(fill = "Dataset Pair", x = "Dataset Pair", y = corr) 

  png(paste0("results/figures/avg/", model, "_model_", corr, "_", scales, "_scales.png"), width=10, height=5, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot model correlations
plot_model(lm_compile, "Linear", "Spearman")
plot_model(lm_compile, "Linear", "Spearman", "free_y")
plot_model(ls_compile, "LASSO", "Spearman")
plot_model(ls_compile, "LASSO", "Spearman", "free_y")

plot_model(lm_compile, "Linear", "Pearson")
plot_model(lm_compile, "Linear", "Pearson", "free_y")
plot_model(ls_compile, "LASSO", "Pearson")
plot_model(ls_compile, "LASSO", "Pearson", "free_y")


############################################################
# Plot feature results
############################################################

# function to plot feature weights
plot_feature <- function(model, label) {

  # subset to model and feature of interest
  toPlot <- feature_df[feature_df$model == model,]
  toPlot <- toPlot[toPlot$label == label,]

  # format dataframe for plotting
  toPlot$Feature <- factor(toPlot$Feature, levels = c("Length", "No. Exons", "GC%", "Median Exp"))
  toPlot$pair <- factor(toPlot$pair, 
                        levels = c("gcsi_ccle", "gcsi_gdsc", "gdsc_ccle"),
                        labels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))

  p <- ggplot(toPlot, aes(x = pair, y = Feature, fill = Weight)) + 
    geom_tile(color = "black") +
    geom_text(aes(label = round(Weight, 3))) +
    scale_fill_gradient2(low = "#9D3737", mid = "white", high = "#3670A0", midpoint = 0) +
    theme_void() + 
    theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust=0.5), 
          axis.text.y = element_text(size = 10, vjust = 0.5, hjust=1), 
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12, angle = 90, vjust = 0.5),
          legend.key.size = unit(0.5, 'cm')) +
    labs(x = "", y = "", fill = "Feature\nWeight", title = label)

  png(paste0("results/figures/features/", model, "_", label, ".png"), width=4.2, height=3, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot feature weights
plot_feature("Linear", "Gene Expression")
plot_feature("Linear", "Isoform Expression")
plot_feature("Linear", "CIRI2")
plot_feature("Linear", "CIRCexplorer2")
plot_feature("Linear", "circRNA_finder")
plot_feature("Linear", "find_circ")

plot_feature("LASSO", "Gene Expression")
plot_feature("LASSO", "Isoform Expression")
plot_feature("LASSO", "CIRI2")
plot_feature("LASSO", "CIRCexplorer2")
plot_feature("LASSO", "circRNA_finder")
plot_feature("LASSO", "find_circ")

############################################################
# FROM HERE: RUN IndivFeatureInfluence.ipynb
############################################################


############################################################
# Define function to load in univariable model results
############################################################

# helper function to format each dataset pair result
load_pair <- function(path, pair, filetype) {

  features <- c("gc", "n_exon", "length")

  # get median
  median <- c("gcsi_ccle" = "gdsc_median", 
              "gcsi_gdsc" = "ccle_median",
              "gdsc_ccle" = "gcsi_median")
  features <- c(median[gsub("/", "", pair)], features)

  # process model datafile
  if (filetype == "model") {

    # create dataframe to store results
    res <- data.frame(matrix(nrow = 0, ncol = 6))
    
    # loop through each feature
    for (feat in features) {
      
      # load in data files
      df <- read.csv(paste0(path, pair, feat, "_lm.csv"))

      # remove brackets from Pearson
      df$Pearson <- gsub("\\[", "", gsub("]", "", df$Pearson)) |> as.numeric()

      # add labels
      df$pair <- gsub("/", "", pair)
      df$feature <- feat
      df$feature[df$feature %in% c("gdsc_median", "ccle_median", "gcsi_median")] <- "median"

      res <- rbind(res, df)
    }
  }

  # process feature datafile
  if (filetype == "features") {

    # create dataframe to store results
    res <- data.frame(matrix(nrow = 0, ncol = 3))

    # loop through each feature
    for (feat in features) {

      # load in data files
      df <- read.csv(paste0(path, pair, feat, "_lm_features.csv"))
      df$pair <- gsub("/", "", pair)

      # rename median
      df$Peak[df$Peak %in% c("gdsc_median", "ccle_median", "gcsi_median")] <- "median"
      colnames(df)[1] <- "Feature"

      res <- rbind(res, df)

    }
  }
  return(res)
}

# function to load in results
load_model <- function(label, filetype) {

  path = paste0("results/data/feature_influence/univariable/", label)

  # process model files
  if (filetype == "model") {
    df <- rbind(
      load_pair(path, "/gcsi_ccle/", "model"),
      load_pair(path, "/gcsi_gdsc/", "model"),
      load_pair(path, "/gdsc_ccle/", "model")
    )
  }
  # process feature files
  if (filetype == "features") {
    df <- rbind(
      load_pair(path, "/gcsi_ccle/", "features"),
      load_pair(path, "/gcsi_gdsc/", "features"),
      load_pair(path, "/gdsc_ccle/", "features")
    )
  }
  df$label <- label
  return(df)
}


############################################################
# Load in univariable results
############################################################

# set up palette for plotting
pal = c("#DACCAB", "#C78B76", "#9D3737")

# set up model results
model_df <- rbind(
  load_model("Gene_Expression", "model"),
  load_model("Isoform_Expression", "model"),
  load_model("CIRI2", "model"),
  load_model("CIRCexplorer2", "model"),
  load_model("circRNA_finder", "model"),
  load_model("find_circ", "model")
)

# format model dataframe
model_df$Fold <- factor(model_df$Fold, levels = c(1, 2, 3, 4, 5))
model_df$pair <- factor(model_df$pair, 
  levels = c("gcsi_ccle", "gcsi_gdsc", "gdsc_ccle"),
  labels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))
model_df$label <- factor(model_df$label, 
  levels = c("Gene_Expression", "Isoform_Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"),
  labels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
model_df$feature <- factor(model_df$feature, 
  levels = c("median", "gc", "n_exon", "length"),
  labels = c("Median Exp", "GC%", "No. Exons", "Length"))

# set up feature results
feature_df <- rbind(
  load_model("Gene_Expression", "features"),
  load_model("Isoform_Expression", "features"),
  load_model("CIRI2", "features"),
  load_model("CIRCexplorer2", "features"),
  load_model("circRNA_finder", "features"),
  load_model("find_circ", "features")
)

# format feature dataframe
feature_df$Feature <- factor(feature_df$Feature, 
  levels = c("median", "gc", "n_exon", "length"),
  labels = c("Median Exp", "GC%", "No. Exons", "Length"))
feature_df$label <- factor(feature_df$label, 
  levels = c("Gene_Expression", "Isoform_Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"),
  labels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))


############################################################
# Format and save model results
############################################################

# function to compile results
compile_res <- function() {
  
  df <- model_df
  df$cat <- paste(df$label, df$pair, df$feature, sep = "-")
  
  # get average and max
  res <- df %>%
    group_by(cat) %>%
    summarise(
      avg_spearman = mean(Spearman, na.rm = TRUE),
      sd_spearman  = sd(Spearman, na.rm = TRUE),
      avg_pearson  = mean(Pearson, na.rm = TRUE),
      sd_pearson   = sd(Pearson, na.rm = TRUE),
      max_spearman = Spearman[which.max(abs(Spearman))],
      max_pearson  = Pearson[which.max(abs(Pearson))],
      .groups = "drop"
    )
  res <- res %>%
    separate(cat, into = c("label", "pair", "feature"), sep = "-", remove = FALSE)
  
  # format results
  res$label <- factor(res$label, levels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
  res$pair <- factor(res$pair, levels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))
  res$feature <- factor(res$feature, levels = c("Median Exp", "GC%", "No. Exons", "Length"))

  res <- res[order(res$label, res$pair, res$feature),] |> as.data.frame()

  # transpose
  rownames(res) <- res$cat
  res <- t(res[,-c(1)]) |> as.data.frame()

  return(res)
}

# compile results
univariable_compile <- compile_res()

# write results
write.csv(univariable_compile, file = "results/data/feature_influence/multivariable/indiv_lm.csv", quote = F, row.names = T)


############################################################
# Plot model results from univariable models (across folds)
############################################################

# function to plot model correlations
plot_model <- function(corr) {

  p <- ggplot(model_df, aes(x = pair, y = .data[[corr]], fill = pair)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    scale_fill_manual(values = pal) +
    facet_nested(~ factor(label) + factor(feature), scales = "free_x") +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm'),
          axis.text.x = element_text(angle = 90, hjust = 0.1)) +
    labs(fill = "Dataset Pair", x = "Dataset Pair") 

  png(paste0("results/figures/univariable_", corr, "_", scales, "_scales.png"), width=12, height=5, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot model correlations
plot_model("Spearman")
plot_model("Pearson")


############################################################
# Plot model results (average across folds)
############################################################

# function to plot model correlations
plot_model <- function(corr, scales = "fixed") {

  toPlot <- t(univariable_compile) |> as.data.frame()

  # format columns
  if (corr == "Spearman") {
    colnames(toPlot)[colnames(toPlot) == "avg_spearman"] <- "corr"
    colnames(toPlot)[colnames(toPlot) == "sd_spearman"] <- "sd"
  }
  if (corr == "Pearson") {
    colnames(toPlot)[colnames(toPlot) == "avg_pearson"] <- "corr"
    colnames(toPlot)[colnames(toPlot) == "sd_pearson"] <- "sd"
  }
  toPlot$corr <- as.numeric(toPlot$corr)
  toPlot$sd <- as.numeric(toPlot$sd)

  toPlot$pair <- factor(toPlot$pair, 
    levels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))
  toPlot$label <- factor(toPlot$label, 
    levels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
  toPlot$feature <- factor(toPlot$feature, levels = c("Median Exp", "GC%", "No. Exons", "Length"))

  p <- ggplot(toPlot, aes(x = feature, y = corr, fill = pair)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=corr-sd, ymax=corr+sd), width=.2, position=position_dodge(.9)) +
    scale_fill_manual(values = pal) +
    facet_wrap(.~label, nrow = 2, scales = scales) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm')) +
    labs(fill = "Dataset Pair", x = "Dataset Pair", y = corr) 

  png(paste0("results/figures/avg/univariable_", corr, "_", scales, "_scales.png"), width=9, height=5, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot model correlations
plot_model("Spearman")
plot_model("Spearman", "free_y")
plot_model("Pearson")
plot_model("Pearson", "free_y")

############################################################
# Plot feature results from univariable models 
############################################################

# function to plot feature weights
plot_feature <- function(scales = "fixed") {

  p <- ggplot(feature_df, aes(x = Feature, y = Weight, fill = pair)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    scale_fill_manual(values = pal, 
                      labels = c("gcsi_ccle" = "gCSI/CCLE",
                                "gcsi_gdsc" = "gCSI/GDSC",
                                "gdsc_ccle" = "GDSC/CCLE")) +
    facet_wrap(.~label, nrow = 1, scales = scales) +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm'),
          axis.text.x = element_text(angle = 90, hjust = 0.1)) +
    labs(fill = "Dataset Pair") 

  png(paste0("results/figures/univariable_features_", scales, "_scales.png"), width=11, height=5, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot feature weights
plot_feature()
plot_feature("free_y")
