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

filename <- paste0("../results/figures/figure7/venndiagram_", analysis, ".png")
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
toPlot$variable <- factor(toPlot$variable, levels = c("count_ribo0", "count_polyA"))
toPlot$value <- log2(toPlot$value + 1)

# plot heatmap of counts
filename <- paste0("../results/figures/figure7/heatmap_", analysis, ".png")
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
  scale_x_discrete(labels = c("rRNA-\ndepleted", "poly(A)\nselected")) +
  guides(fill = guide_colourbar(barwidth = 0.2, barheight = 2, title = "log2(Counts)")) +
  labs(x = "\nPipeline", y = "Tumour Sample")
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

filename <- paste0("../results/figures/figure7/prop_unique_", analysis, ".png")
print(paste("Saving donuts to", filename))
png(filename, width=12, height=3, units='in', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = "right")
dev.off()

print("done")