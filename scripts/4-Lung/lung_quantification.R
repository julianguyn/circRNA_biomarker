# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(viridis)
    library(reshape2)
    library(ggvenn)
})

options(stringsAsFactors = FALSE)
source("utils/palettes.R")

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
# Venn diagram of unique transcripts per protocol
############################################################

# helper function to plot venn diagram
plot_venn_protocol <- function(
  ciri, circ, cfnd, fcrc, label
) {
  transcripts <- list(
    CIRI2 = colnames(ciri),
    CIRCexplorer2 = colnames(circ),
    circRNA_finder = colnames(cfnd),
    find_circ = colnames(fcrc)
  )
  p <- ggvenn(transcripts,fill_color = pipeline_pal, stroke_size = 0.5, set_name_size = 4) +
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    labs(title = paste0(label, "\n"))
  return(p)
}

p1 <- plot_venn_protocol(ciri_polyA, circ_polyA, cfnd_polyA, fcrc_polyA, "Poly(A)-Selection")
p2 <- plot_venn_protocol(ciri_ribo0, circ_ribo0, cfnd_ribo0, fcrc_ribo0, "rRNA-Depletion")

filename <- paste0("../results/figures/figure4/venndiagram_", analysis, ".png")
png(filename, width=250, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = FALSE)
dev.off()

############################################################
# Plot heatmap of counts
############################################################

# helper function
combine_pipeline <- function(polyA, ribo0, pipeline) {
  df <- data.frame(
    samples = paste0("tumour", c(1:51)),
    count_polyA = rowSums(polyA),
    count_ribo0 = rowSums(ribo0)
  )
  df$pipeline <- pipeline
  rownames(df) <- NULL
  return(df)
}

combined_df <- rbind(
  combine_pipeline(ciri_polyA, ciri_ribo0, "CIRI2"),
  combine_pipeline(circ_polyA, circ_ribo0, "CIRCexplorer2"),
  combine_pipeline(cfnd_polyA, cfnd_ribo0, "circRNA_finder"),
  combine_pipeline(fcrc_polyA, fcrc_ribo0, "find_circ")
)

# format dataframe for plotting
toPlot <- reshape2::melt(combined_df)
toPlot$variable <- factor(toPlot$variable, levels = c("count_polyA", "count_ribo0"))
toPlot$value <- log2(toPlot$value + 1)

# plot heatmap of counts
filename <- paste0("../results/figures/figure4/heatmap_", analysis, ".png")
print(paste("Saving heatmap to", filename))
png(filename, width=250, height=200, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = variable, y = samples, fill = value)) + 
  geom_tile() + 
  facet_grid(~factor(pipeline, levels=c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))) +
  theme_void() +
  theme(text = element_text(size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17, angle = 90),
        strip.text.x = element_text(size = 17),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  scale_fill_viridis(limits=c(0, 13), option="mako", direction = -1) +
  scale_x_discrete(labels = c("poly(A)\nselected", "rRNA-\ndepleted")) +
  guides(fill = guide_colourbar(barwidth = 0.2, barheight = 2, title = "log2(Counts)")) +
  labs(x = "\nProtocol", y = "Tumour Sample")
dev.off()


############################################################
# Plot proportion of unique transcripts per pipeline
############################################################

# function to compute proportion of unique transcripts per pipeline
plot_unique_transcripts <- function(polyA_df, ribo0_df, title) {

  # common transcripts
  transcripts <- intersect(colnames(polyA_df), colnames(ribo0_df))

  # count of unique transcript detection per protocol
  count <- data.frame(
    protocol = c("polyA", "ribo0", "both"),
    count = c(
      length(colnames(polyA_df)) - length(transcripts),
      length(colnames(ribo0_df)) - length(transcripts),
      length(transcripts)
    )
  )
  # formating for plot
  count$fraction <- count$count / sum(count$count)
  count$ymax <- cumsum(count$fraction)
  count$ymin <- c(0, head(count$ymax, n = -1))

  # make label for plot
  pt <- round(count$fraction[count$protocol == "both"] * 100, 2)
  label <- paste0(length(transcripts), "\n(", pt, "%)")

  p <- ggplot(count, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = protocol)) +
    geom_rect(color = "black") +
    coord_polar(theta="y") +
    scale_fill_manual(
      values = protocol_pal,
      labels = c(
        "Both protocols",
        "poly(A)-selected only     ",
        "rRNA-depleted only"
      )) + 
    xlim(c(2, 4)) + 
    annotate("text", x = 2, y = 0.85, label = label, size = 5) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      text = element_text(size = 13),
      legend.key.size = unit(1, 'cm'), 
      legend.position = "bottom", 
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.text.x = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    ) +
    labs(title = title, fill = "Transcript detection")
  
  return(p)
}

# plot donut plots
p1 <- plot_unique_transcripts(ciri_polyA, ciri_ribo0, "CIRI2")
p2 <- plot_unique_transcripts(circ_polyA, circ_ribo0, "CIRCexplorer2")
p3 <- plot_unique_transcripts(cfnd_polyA, cfnd_ribo0, "circRNA_finder")
p4 <- plot_unique_transcripts(fcrc_polyA, fcrc_ribo0, "find_circ")

filename <- paste0("../results/figures/figure4/prop_unique_", analysis, ".png")
print(paste("Saving donuts to", filename))
png(filename, width=12, height=3, units='in', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = "right")
dev.off()

# -------------- Pick one plot:

############################################################
# Plot difference in expression between matched transcripts
############################################################

diff_common_transcript <- function(polyA_df, ribo0_df, pipeline) {

  # keep common transcripts
  transcripts <- intersect(colnames(polyA_df), colnames(ribo0_df))
  polyA_df <- polyA_df[, transcripts]
  ribo0_df <- ribo0_df[, transcripts]

  # get difference of log2 expression
  polyA_df <- log2(polyA_df + 0.001)
  ribo0_df <- log2(ribo0_df + 0.001)

  diff_df <- ribo0_df - polyA_df

  # remove those where expression was 0 in both datasets
  keep <- !(polyA_df == 0 & ribo0_df == 0)
  diff <- diff_df[keep]

  # make results
  df <- data.frame(pipeline = pipeline, diff = diff)
  return(df)
}

toPlot <- rbind(
  diff_common_transcript(ciri_polyA, ciri_ribo0, "CIRI2"),
  diff_common_transcript(circ_polyA, circ_ribo0, "CIRCexplorer2"),
  diff_common_transcript(cfnd_polyA, cfnd_ribo0, "circRNA_finder"),
  diff_common_transcript(fcrc_polyA, fcrc_ribo0, "find_circ")
)
toPlot$pipeline <- factor(toPlot$pipeline, levels=c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

# plot density plots of differences
filename <- paste0("../results/figures/figure4/density_", analysis, ".png")
print(paste("Saving density to", filename))
png(filename, width=7, height=3, units='in', res = 600, pointsize=80)
ggplot(toPlot, aes(x = diff, fill = pipeline)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = pipeline_pal) +
  theme_classic() +
  theme(legend.key.size = unit(0.6, 'cm')) +
  labs(x = "\nDifference in log2(expression)", y = "Density", fill = "Pipeline")
dev.off()

# plot violin plots of differences
filename <- paste0("../results/figures/figure4/violin_", analysis, ".png")
print(paste("Saving violin to", filename))
png(filename, width=7, height=3, units='in', res = 600, pointsize=80)
ggplot(toPlot, aes(x = diff, y = pipeline, fill = pipeline)) +
  geom_violin() +
  scale_fill_manual(values = pipeline_pal) +
  theme_classic() +
  theme(legend.key.size = unit(0.6, 'cm')) +
  labs(x = "\nDifference in log2(expression)", y = NULL, fill = "Pipeline")
dev.off()


############################################################
# Plot heatmap of differences
############################################################

range <- max(max(toPlot$diff, abs(min(toPlot$diff)))) |> round() + 1

plot_common_transcript <- function(polyA_df, ribo0_df, pipeline) {

  # keep common transcripts
  transcripts <- intersect(colnames(polyA_df), colnames(ribo0_df))
  polyA_df <- polyA_df[, transcripts]
  ribo0_df <- ribo0_df[, transcripts]

  # get difference of log2 expression
  polyA_df <- log2(polyA_df + 0.001)
  ribo0_df <- log2(ribo0_df + 0.001)

  diff_df <- ribo0_df - polyA_df

  # remove those where expression was 0 in both datasets
  keep <- !(polyA_df == 0 & ribo0_df == 0)
  diff <- diff_df[keep]

  # make results
  df <- data.frame(pipeline = pipeline, diff = diff)
  df <- df[order(df$diff, decreasing = TRUE), ]
  df$rank <- 1:nrow(df)
  
  p <- ggplot(df, aes(y = 1, x = rank, fill = diff)) +
    geom_tile(linewidth = 0) +
    scale_fill_gradient2(
      "Difference\nin log2(exp)",
      low = "#BC4749",
      high = "#689CB0",
      mid = "#C2BBC9",
      limits = c(-range, range)
    ) +
    theme_void() +
    theme(
        axis.title.y = element_text(size=9, hjust=1),
        legend.title = element_text(size=8),
        axis.ticks.y = element_line(color = "gray", linewidth = 0.3)
    ) + 
    labs(y = pipeline)
  return(p)
}

# plot donut plots
p1 <- plot_common_transcript(ciri_polyA, ciri_ribo0, "c1")
p2 <- plot_common_transcript(circ_polyA, circ_ribo0, "c2")
p3 <- plot_common_transcript(cfnd_polyA, cfnd_ribo0, "c3")
p4 <- plot_common_transcript(fcrc_polyA, fcrc_ribo0, "c4")

filename <- paste0("../results/figures/figure4/tile_", analysis, ".png")
print(paste("Saving tiles to", filename))
png(filename, width=6.5, height=2.5, units='in', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, nrow = 4, common.legend = TRUE, legend = "right")
dev.off()

print("done")