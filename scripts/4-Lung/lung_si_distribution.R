# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(ggvenn)
    library(reshape2)
    library(matrixStats)
})

options(stringsAsFactors = FALSE)
set.seed(123)

source("utils/palettes.R")
source("utils/compute_stability_index.R")

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("circ", "GE")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

############################################################
# Load in and prepare metadata
############################################################

# load in metadata
polyA_meta <- read.table("../data/rnaseq_meta/lung_polyA.tsv", header = TRUE)
ribo0_meta <- read.csv("../data/rnaseq_meta/lung_ribozero.csv")

# keep common 51 (and order)
polyA_meta <- polyA_meta[match(ribo0_meta$TB_id, polyA_meta$TB_id), ]


############################################################
# Load in counts matrices and keep common samples
############################################################

# helper function to load in count matrices
load_lung <- function(filename, analysis) {
  # load data
  indir <- paste0("../data/processed_lung/", analysis, "/")
  df <- fread(paste0(indir, filename), data.table = FALSE)
  # get protocol
  protocol <- sub(".*_", "", sub("_counts.tsv", "", filename))
  # match samples to metadata
  if (protocol == "polyA") {
    df <- df[match(polyA_meta$sample, df$sample), ]
  } else if (protocol == "ribo0") {
    df <- df[match(ribo0_meta$helab_id, df$sample), ]
  } else {
    print(paste("Unable to determine protocol for", filename))
  }
  rownames(df) <- paste0("tumour", c(1:51))
  df$sample <- NULL

  # filter to keep only transcripts in at least 6 samples
  df <- df[,-which(colSums(df == 0) > 45)]
  print(dim(df))
  return(df)
}

indir <- paste0("../data/processed_lung/", analysis, "/")
print(paste("Loading in files from:", indir))

# load in count matrices
ciri_polyA <- load_lung("ciri_polyA_counts.tsv", analysis)
ciri_ribo0 <- load_lung("ciri_ribo0_counts.tsv", analysis)
circ_polyA <- load_lung("circ_polyA_counts.tsv", analysis)
circ_ribo0 <- load_lung("circ_ribo0_counts.tsv", analysis)
cfnd_polyA <- load_lung("cfnd_polyA_counts.tsv", analysis)
cfnd_ribo0 <- load_lung("cfnd_ribo0_counts.tsv", analysis)
fcrc_polyA <- load_lung("fcrc_polyA_counts.tsv", analysis)
fcrc_ribo0 <- load_lung("fcrc_ribo0_counts.tsv", analysis)


############################################################
# Compute pairwise Spearman corr from datasets
############################################################

print("Starting polyA stability indices")
polyA_stability <- rbind(
    compute_stability_index(ciri_polyA, circ_polyA, "ciri", "circ", "polyA"),
    compute_stability_index(ciri_polyA, cfnd_polyA, "ciri", "cfnd", "polyA"),
    compute_stability_index(ciri_polyA, fcrc_polyA, "ciri", "fcrc", "polyA"),
    compute_stability_index(circ_polyA, cfnd_polyA, "circ", "cfnd", "polyA"),
    compute_stability_index(circ_polyA, fcrc_polyA, "circ", "fcrc", "polyA"),
    compute_stability_index(cfnd_polyA, fcrc_polyA, "cfnd", "fcrc", "polyA")
)

print("Starting ribo0 stability indices")
ribo0_stability <- rbind(
    compute_stability_index(ciri_ribo0, circ_ribo0, "ciri", "circ", "ribo0"),
    compute_stability_index(ciri_ribo0, cfnd_ribo0, "ciri", "cfnd", "ribo0"),
    compute_stability_index(ciri_ribo0, fcrc_ribo0, "ciri", "fcrc", "ribo0"),
    compute_stability_index(circ_ribo0, cfnd_ribo0, "circ", "cfnd", "ribo0"),
    compute_stability_index(circ_ribo0, fcrc_ribo0, "circ", "fcrc", "ribo0"),
    compute_stability_index(cfnd_ribo0, fcrc_ribo0, "cfnd", "fcrc", "ribo0")
)

print("Starting polyA random stability indices")
polyA_stability_random <- rbind(
    compute_stability_index(ciri_polyA, circ_polyA, "ciri", "circ", "polyA", random = TRUE, iter = 100),
    compute_stability_index(ciri_polyA, cfnd_polyA, "ciri", "cfnd", "polyA", random = TRUE, iter = 100),
    compute_stability_index(ciri_polyA, fcrc_polyA, "ciri", "fcrc", "polyA", random = TRUE, iter = 100),
    compute_stability_index(circ_polyA, cfnd_polyA, "circ", "cfnd", "polyA", random = TRUE, iter = 100),
    compute_stability_index(circ_polyA, fcrc_polyA, "circ", "fcrc", "polyA", random = TRUE, iter = 100),
    compute_stability_index(cfnd_polyA, fcrc_polyA, "cfnd", "fcrc", "polyA", random = TRUE, iter = 100)
)

print("Starting ribo0 random stability indices")
ribo0_stability_random <- rbind(
    compute_stability_index(ciri_ribo0, circ_ribo0, "ciri", "circ", "ribo0", random = TRUE, iter = 100),
    compute_stability_index(ciri_ribo0, cfnd_ribo0, "ciri", "cfnd", "ribo0", random = TRUE, iter = 100),
    compute_stability_index(ciri_ribo0, fcrc_ribo0, "ciri", "fcrc", "ribo0", random = TRUE, iter = 100),
    compute_stability_index(circ_ribo0, cfnd_ribo0, "circ", "cfnd", "ribo0", random = TRUE, iter = 100),
    compute_stability_index(circ_ribo0, fcrc_ribo0, "circ", "fcrc", "ribo0", random = TRUE, iter = 100),
    compute_stability_index(cfnd_ribo0, fcrc_ribo0, "cfnd", "fcrc", "ribo0", random = TRUE, iter = 100)
)

# save checkpoint
filename <- paste0("../results/data/temp/lung_stability_", analysis, ".RData")
save(polyA_stability, ribo0_stability, polyA_stability_random, ribo0_stability_random, file = filename)


############################################################
# Format stability index matrices for plotting
############################################################

# function to format stability dataframes
format_df <- function(df, label, random = "NonRandom") {
    toPlot <- reshape2::melt(df)
    colnames(toPlot) <- c("Pipeline", "Stability")
    toPlot$label <- label
    toPlot$random <- random
    return(toPlot)
}

polyA_stability <- format_df(polyA_stability, "polyA", "NonRandom")
ribo0_stability <- format_df(ribo0_stability, "ribo0", "NonRandom")

polyA_stability_random <- format_df(polyA_stability_random, "polyA", "Random")
ribo0_stability_random <- format_df(ribo0_stability_random, "ribo0", "Random")


############################################################
# Plot stability index distribution 
############################################################

# merge nonrandom results for plotting
toPlot <- rbind(polyA_stability, ribo0_stability, polyA_stability_random, ribo0_stability_random)
toPlot$label <- factor(toPlot$label, levels = c("polyA", "ribo0"), labels = c("PolyA", "Ribo0"))

png("../results/figures/figure7/stability_lung.png", width=100, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = Stability, fill = label)) + 
    geom_boxplot(data = toPlot, aes(alpha = random)) + 
    scale_fill_manual(values = c("#4CC5AB", "#392C57", "grey")) +
    scale_alpha_manual(values = c(1, 0.2)) +
    facet_grid(factor(Pipeline)~.) +
    theme_classic() + labs(x = "", fill = "", y = "Stability Index", alpha = "Randomization") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()
