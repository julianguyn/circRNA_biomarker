# load libraries
suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(tidyverse)
    library(ggh4x)
    library(RColorBrewer)
    library(PharmacoGx)
})

set.seed(200)  


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
    lm$l1_ratio <- lm$alpha <- lm$max_iter <- NA
    ls <- read.csv(paste0(path, pair, "lasso.csv"))
    ls$l1_ratio <- NA
    en <- read.csv(paste0(path, pair, "en.csv"))
    en$Model <- "ElasticNet"

    # merge data files
    df <- rbind(lm, ls, en)
    df$pair <- gsub("/", "", pair)
  }

  # process feature datafile
  if (filetype == "features") {

    # load in data files
    lm <- read.csv(paste0(path, pair, "lm_features.csv"))
    ls <- read.csv(paste0(path, pair, "lasso_features.csv"))
    en <- read.csv(paste0(path, pair, "en_features.csv"))

    # add in labels
    lm$model <- "Linear"
    ls$model <- "LASSO"
    en$model <- "ElasticNet"

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

  path <- paste0("results/data/feature_influence/multivariable/", label)

  df <- switch(
    filetype,
    model = rbind(
      load_pair(path, "/gcsi_ccle/", "model"),
      load_pair(path, "/gcsi_gdsc/", "model"),
      load_pair(path, "/gdsc_ccle/", "model")
    ),
    features = rbind(
      load_pair(path, "/gcsi_ccle/", "features"),
      load_pair(path, "/gcsi_gdsc/", "features"),
      load_pair(path, "/gdsc_ccle/", "features")
    )
  )
  df$label <- label
  return(df)
}


############################################################
# Load in mulivariable results
############################################################

# set up palette for plotting
pal = c("#DACCAB", "#C78B76", "#9D3737", "gray")

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
      ci_spearman  = sd(Spearman, na.rm = TRUE)/ sqrt(5) * 1.96,
      avg_pearson  = mean(Pearson, na.rm = TRUE),
      ci_pearson   = sd(Pearson, na.rm = TRUE)/ sqrt(5) * 1.96,
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
en_compile <- compile_res("ElasticNet")

# write results
write.csv(lm_compile, file = "results/data/feature_influence/multivariable/lm.csv", quote = F, row.names = T)
write.csv(ls_compile, file = "results/data/feature_influence/multivariable/LASSO.csv", quote = F, row.names = T)
write.csv(en_compile, file = "results/data/feature_influence/multivariable/en.csv", quote = F, row.names = T)


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
plot_model("ElasticNet", "Spearman")
plot_model("ElasticNet", "Spearman", "free_y")

plot_model("Linear", "Pearson")
plot_model("Linear", "Pearson", "free_y")
plot_model("LASSO", "Pearson")
plot_model("LASSO", "Pearson", "free_y")
plot_model("ElasticNet", "Pearson")
plot_model("ElasticNet", "Pearson", "free_y")


############################################################
# Plot model results (average across folds)
############################################################

# function to plot model correlations
plot_model <- function(df, model, corr, scales = "fixed") {

  toPlot <- t(df) |> as.data.frame()

  # format columns
  if (corr == "Spearman") {
    colnames(toPlot)[colnames(toPlot) == "avg_spearman"] <- "corr"
    colnames(toPlot)[colnames(toPlot) == "ci_spearman"] <- "ci"
  }
  if (corr == "Pearson") {
    colnames(toPlot)[colnames(toPlot) == "avg_pearson"] <- "corr"
    colnames(toPlot)[colnames(toPlot) == "ci_pearson"] <- "ci"
  }
  toPlot$corr <- as.numeric(toPlot$corr)
  toPlot$ci <- as.numeric(toPlot$ci)

  toPlot$pair <- gsub(".*-", "", rownames(toPlot))
  toPlot$label <- gsub("-.*", "", rownames(toPlot))

  toPlot$label <- factor(toPlot$label, 
    levels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

  # set fill colours for CI past zero
  toPlot <- toPlot %>%
    mutate(crosses_zero = (corr - ci <= 0 & corr + ci >= 0))
  toPlot <- toPlot %>%
    mutate(fill = ifelse(crosses_zero, "Non-Significant", pair))
  toPlot$fill[is.na(toPlot$fill)] <- "Non-Significant"
  toPlot$fill <- factor(toPlot$fill, 
    levels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE", "Non-Significant"))


  p <- ggplot(toPlot, aes(x = pair, y = corr, fill = fill)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=corr-ci, ymax=corr+ci), width=.2, position=position_dodge(.9)) +
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

  if (scales == "fixed") return(toPlot)
}

# plot model correlations
lm_s <- plot_model(lm_compile, "Linear", "Spearman")
plot_model(lm_compile, "Linear", "Spearman", "free_y")
ls_s <- plot_model(ls_compile, "LASSO", "Spearman")
plot_model(ls_compile, "LASSO", "Spearman", "free_y")
en_s <- plot_model(en_compile, "ElasticNet", "Spearman")
plot_model(en_compile, "ElasticNet", "Spearman", "free_y")

lm_p <- plot_model(lm_compile, "Linear", "Pearson")
plot_model(lm_compile, "Linear", "Pearson", "free_y")
ls_p <- plot_model(ls_compile, "LASSO", "Pearson")
plot_model(ls_compile, "LASSO", "Pearson", "free_y")
en_p <- plot_model(en_compile, "ElasticNet", "Pearson")
plot_model(en_compile, "ElasticNet", "Pearson", "free_y")

############################################################
# Plot all model results
############################################################

toPlot <- rbind(lm_s, ls_s, en_s)
toPlot$model <- rep(c("Linear", "LASSO", "ElasticNet"), each = 18)

p <- ggplot(toPlot, aes(x = pair, y = corr, fill = model)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=corr-ci, ymax=corr+ci), width=.2, position=position_dodge(.9)) +
    scale_fill_manual(values = pal) +
    facet_wrap(.~label, nrow = 1, scales = "fixed") +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm'),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(fill = "Dataset Pair", x = "Dataset Pair", y = corr)

png("results/figures/all_models.png", width=10, height=4, units='in', res = 600, pointsize=80)
p
dev.off()

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

plot_feature("ElasticNet", "Gene Expression")
plot_feature("ElasticNet", "Isoform Expression")
plot_feature("ElasticNet", "CIRI2")
plot_feature("ElasticNet", "CIRCexplorer2")
plot_feature("ElasticNet", "circRNA_finder")
plot_feature("ElasticNet", "find_circ")


############################################################
# Feature distribution by SI quantiles
############################################################

# load in stability indices
gene <- read.csv("temp/gene_stability.csv")
isof <- read.csv("temp/transcript_stability.csv")
ciri <- read.csv("temp/ciri_stability.csv")
circ <- read.csv("temp/circ_stability.csv")
cfnd <- read.csv("temp/cfnd_stability.csv")
fcrc <- read.csv("temp/fcrc_stability.csv")

# helper function to create feature df per pair
helper_SI_quantiles <- function(df, pair) {

  median <- c("gcsi_ccle" = "gdsc_median",
              "gcsi_gdsc" = "ccle_median",
              "gdsc_ccle" = "gcsi_median")

  # columns to keep
  keep <- c(paste0(pair, "_spearman"), "gc", "n_exon", "length", median[pair])
  subset <- df[,colnames(df) %in% keep]
  colnames(subset)[colnames(subset) == median[pair]] <- "median"

  # get quantile labels
  subset$qt <- cut(subset[[paste0(pair, "_spearman")]],
      breaks = quantile(subset[[paste0(pair, "_spearman")]], probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = 1:4)
  subset[[paste0(pair, "_spearman")]] <- NULL

  # melt dataframe
  df <- reshape2::melt(subset, id.vars = "qt")
  df$pair <- pair

  return(df)
}

# function to stratify and plot features per SI quantile
plot_SI_quantiles <- function(df, label) {

  toPlot <- rbind(
    helper_SI_quantiles(df, "gcsi_ccle"),
    helper_SI_quantiles(df, "gcsi_gdsc"),
    helper_SI_quantiles(df, "gdsc_ccle")
  )

  toPlot$variable <- factor(toPlot$variable, 
    levels = c("median", "gc", "n_exon", "length"),
    labels = c("Median Expression", "GC%", "Number of Exons", "Length"))
  toPlot$pair <- factor(toPlot$pair, 
    levels = c("gcsi_ccle", "gcsi_gdsc", "gdsc_ccle"),
    labels = c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE"))
  
  # box plots
  p <- ggplot(toPlot, aes(x = pair, y = value, fill = qt)) + 
    geom_boxplot() + 
    facet_wrap(.~variable, scales = "free_y") + 
    scale_fill_manual(values = c("#333d29", "#414833", "#656d4a", "#a4ac86")) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm')) +
    labs(y = "Feature Value", x = "", fill = "Quantile", title = label)
  
  # save plot
  filepath = paste0("results/figures/distribution/si_", label, "_feature_distribution.png")
  png(filepath, width=180, height=150, units='mm', res = 600, pointsize=80)
  print({p})
  dev.off()

  # return dataframe with median expression only
  toPlot <- toPlot[toPlot$variable == "Median Expression",]
  toPlot$label <- label
  return(toPlot)
}

# compile dataframes (and plot)
toPlot <- rbind(
  plot_SI_quantiles(gene, "Gene Expression"),
  plot_SI_quantiles(isof, "Isoform Expression"),
  plot_SI_quantiles(ciri, "CIRI2"),
  plot_SI_quantiles(circ, "CIRCexplorer2"),
  plot_SI_quantiles(cfnd, "circRNA_finder"),
  plot_SI_quantiles(fcrc, "find_circ")
)
toPlot$label <- factor(toPlot$label, 
    levels = c("Gene Expression", "Isoform Expression", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
  

# plot median expression 
toPlot <- toPlot[!is.na(toPlot$qt),]
png("results/figures/distribution/median_expression.png", width=13.5, height=5, units='in', res = 600, pointsize=80)
ggplot(toPlot, aes(x = pair, y = value, fill = qt)) + 
  geom_boxplot() + 
  facet_wrap(.~label, scales = "free_y", nrow = 1) + 
  scale_fill_manual(values = c("#333d29", "#414833", "#656d4a", "#a4ac86")) +
  theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(angle = 90, hjust = 0.1)) +
  labs(y = "Median Expression", x = "", fill = "Quantile")
dev.off()
