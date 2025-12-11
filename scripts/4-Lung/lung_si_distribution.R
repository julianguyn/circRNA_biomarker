# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(ggvenn)
    library(reshape2)
    library(matrixStats)
    library(stats)
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

  # filter to keep only transcripts in at least 6 samples **************
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
# Plot stability index distribution
############################################################

# merge nonrandom results for plotting
toPlot <- rbind(polyA_stability, ribo0_stability, polyA_stability_random, ribo0_stability_random)
toPlot$dataset <- factor(toPlot$dataset, levels = c("polyA", "ribo0"))
toPlot$label <- factor(
  toPlot$label,
  levels = c("ciri/circ", "ciri/cfnd", "ciri/fcrc", "circ/cfnd", "circ/fcrc", "cfnd/fcrc")
)

filename <- paste0("../results/figures/figure4/stability_lung_", analysis, ".png")
png(filename, width=4, height=10, units='in', res = 600, pointsize=80)
ggplot(toPlot, aes(x = dataset, y = stability, fill = dataset)) +
    geom_boxplot(data = toPlot, aes(alpha = random)) +
    scale_fill_manual(values = protocol_pal, labels = c("poly(A)-selection", "rRNA-depletion")) +
    scale_alpha_manual(values = c(1, 0.2)) +
    facet_grid(factor(label)~.) +
    scale_x_discrete(labels = c("poly(A)\nselected", "rRNA-\ndepleted")) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "", fill = "", y = "Stability Index", alpha = "Randomization") +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()

############################################################
# Wilcoxon rank sum test (+ numbers)
############################################################

for (pair in unique(toPlot$label)) {
  print(paste("Working on:", pair))
  subset <- toPlot[toPlot$label == pair, ]

  p_ran <- subset[which(subset$dataset == "polyA" & subset$random == "Random"), ]
  p_non <- subset[which(subset$dataset == "polyA" & subset$random == "NonRandom"), ]
  r_ran <- subset[which(subset$dataset == "ribo0" & subset$random == "Random"), ]
  r_non <- subset[which(subset$dataset == "ribo0" & subset$random == "NonRandom"), ]

  # print numbers
  print(paste("---", "PolyA:", nrow(p_non)))
  print(paste("---", "Ribo0:", nrow(r_non)))

  if (nrow(p_non) == 0 || nrow(r_non) == 0) {
    message("")
  } else {
    # all nonrandom vs random
    p_w <- wilcox.test(p_non$stability, p_ran$stability, alternative = "greater", exact = FALSE)
    r_w <- wilcox.test(r_non$stability, r_ran$stability, alternative = "greater", exact = FALSE)
    p <- ifelse(p_w$p.value < 0.05, "*", " ")
    r <- ifelse(r_w$p.value < 0.05, "*", " ")
    print("Nonrandom > random:")
    print(paste("---", "polya:", p_w$p.value, p))
    print(paste("---", "ribo0:", r_w$p.value, r))

    # all polyA vs ribo0
    w <- wilcox.test(r_non$stability, p_non$stability, alternative = "greater", exact = FALSE)
    p <- ifelse(w$p.value < 0.05, "*", " ")
    print("Ribo0 > polyA:")
    print(paste("---", w$p.value, p))
    message("")
  }
}