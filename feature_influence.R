# load libraries
suppressPackageStartupMessages({
    library(caret)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(tidyverse)
    library(ggh4x)
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
# Plot model results from linear model
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
# Plot feature results from linear model
############################################################

# function to plot feature weights
plot_feature <- function(model, scales = "fixed") {

  toPlot <- feature_df[feature_df$model == model,]

  p <- ggplot(toPlot, aes(x = Feature, y = Weight, fill = pair)) + 
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

  png(paste0("results/figures/", model, "_features_", scales, "_scales.png"), width=11, height=5, units='in', res = 600, pointsize=80)
  print({p})
  dev.off()
}

# plot feature weights
plot_feature("Linear")
plot_feature("Linear", "free_y")
plot_feature("LASSO")
plot_feature("LASSO", "free_y")


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
# Plot model results from linear model
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
# Plot feature results from linear model
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
