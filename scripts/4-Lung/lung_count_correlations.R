# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(viridis)
    library(reshape2)
    library(patchwork)
})

options(stringsAsFactors = FALSE)
source("scripts/utils/palettes.R")

analysis <- "circ"

############################################################
# Load in and prepare metadata
############################################################

# load in metadata
polyA_meta <- read.table("data/rnaseq_meta/lung_polyA.tsv", header = TRUE)
ribo0_meta <- read.csv("data/rnaseq_meta/lung_ribozero.csv")

# keep common 51 (and order)
polyA_meta <- polyA_meta[match(ribo0_meta$TB_id, polyA_meta$TB_id), ]

polyA_meta$tumourID <- paste0("tumour", c(1:51))
ribo0_meta$tumourID <- paste0("tumour", c(1:51))

############################################################
# Load in rRNA contamination reports
############################################################

load_reports <- function(protocol) {
    path <- paste0("results/data/", protocol, "_rRNA_contamination_summary.csv")
    df <- read.csv(path)
    if (protocol == "polyA") {
        df$Sample <- sub("^([^_]+_[^_]+_[^_]+)_.*", "\\1", sub(".*/", "", df$Sample))
        df <- df[df$Sample %in% polyA_meta$sample,]
        df$tumourID <- polyA_meta$tumourID[match(df$Sample, polyA_meta$sample)]
    } else if (protocol == "ribo0") {
        df$Sample <- sub("_.*", "", sub(".*/", "", df$Sample))
        df$tumourID <- ribo0_meta$tumourID[match(df$Sample, ribo0_meta$helab_id)]
    }
    return(df)
}

polyA_rRNA <- load_reports("polyA")
ribo0_rRNA <- load_reports("ribo0")

############################################################
# Load in sequencing depth
############################################################

# load in seqdepth
polyA_lib <- read.table("data/raw_lung/polyA.tsv", header = TRUE)
polyA_lib <- polyA_lib[polyA_lib$sample %in% polyA_meta$sample,]
polyA_lib$tumourID <- polyA_meta$tumourID[match(polyA_lib$sample, polyA_meta$sample)]

ribo0_lib <- read.table("data/raw_lung/ribo0.tsv", header = TRUE)
ribo0_lib$tumourID <- ribo0_meta$tumourID[match(ribo0_lib$sample, ribo0_meta$helab_id)]

############################################################
# Load in counts matrices and keep common samples
############################################################

# helper function to load in count matrices
load_lung <- function(filename, analysis) {
  # load data
  indir <- paste0("data/processed_lung/", analysis, "/")
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
# Format counts for plotting
############################################################

# helper function
combine_pipeline <- function(ciri, circ, cfnd, fcrc) {
  df <- data.frame(
    sample = paste0("tumour", c(1:51)),
    CIRI2 = rowSums(ciri),
    CIRCexplorer2 = rowSums(circ),
    circRNA_finder = rowSums(cfnd),
    find_circ = rowSums(fcrc)
  )
  rownames(df) <- NULL

  toPlot <- reshape2::melt(df)
  colnames(toPlot) <- c("sample", "pipeline", "count")
  toPlot$pipeline <- factor(toPlot$pipeline, levels = names(pipeline_pal)[1:4])
  toPlot$count <- log2(toPlot$count + 1)
  return(toPlot)
}

polyA_df <- combine_pipeline(ciri_polyA, circ_polyA, cfnd_polyA, fcrc_polyA)
ribo0_df <- combine_pipeline(ciri_ribo0, circ_ribo0, cfnd_ribo0, fcrc_ribo0)


############################################################
# Correlations (CPM and rRNA contamination & seqdepth)
############################################################

get_correlations <- function(toPlot, rRNA_df, lib_df, label) {

    res <- data.frame(matrix(nrow=0, ncol=3))

    for (pipeline in unique(toPlot$pipeline)) {
        subset_df <- toPlot[toPlot$pipeline == pipeline,]
        subset_df$rRNA <- rRNA_df$rRNA_Contamination_Pct[match(subset_df$sample, rRNA_df$tumourID)]
        subset_df$lib <- lib_df$avg_counts[match(subset_df$sample, lib_df$tumourID)]

        df <- data.frame(
            pipeline = pipeline,
            rRNA = suppressWarnings(cor.test(subset_df$count, subset_df$rRNA, method = "spearman")$estimate),
            lib = suppressWarnings(cor.test(subset_df$count, subset_df$lib, method = "spearman")$estimate)
        )
        
        res <- rbind(res, df)
    }
    res$protocol <- label
    return(res)
}

toPlot <- rbind(
    get_correlations(polyA_df, polyA_rRNA, polyA_lib, "polyA"),
    get_correlations(ribo0_df, ribo0_rRNA, ribo0_lib, "ribo0")
)
toPlot$pipeline <- factor(toPlot$pipeline, levels = rev(names(pipeline_pal)[1:4]))


plot_correlation <- function(feature, label) {
    p <- ggplot(toPlot, aes(y = pipeline, x = protocol, fill = .data[[feature]])) +
        geom_tile(color = "gray") +
        geom_text(aes(label = round(.data[[feature]], 2)), size = 4) +
        scale_fill_gradient2("Spearman\nCorrelation", high = "blue", low = "red", mid = "white", limits = c(-1, 1)) +
        scale_x_discrete(labels = c("poly(A)\nselection", "rRNA-\ndepletion")) +
        theme_minimal() +
        theme(
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 11)
        ) +
        labs(title = label)
    return(p)
}

p1 <- plot_correlation("rRNA", "rRNA Contamination")
p2 <- plot_correlation("lib", "Sequencing Depth")

p <- p1 + p2 + plot_layout(guides = "collect")
filename <- "results/figures/suppfig5/correlations.png"
ggsave(filename, p, width = 6, height = 3)


############################################################
# Plot rRNA contamination and counts
############################################################

plot_rRNA_contamination <- function(toPlot, rRNA_df, label) {

  summary_df <- rRNA_df %>%
    group_by(tumourID) %>%
    summarise(mean_pct = mean(rRNA_Contamination_Pct), .groups = "drop")
  summary_df <- summary_df[order(summary_df$mean_pct),]

  summary_df$tumourID <- factor(summary_df$tumourID, levels = summary_df$tumourID)
  rRNA_df$tumourID <- factor(rRNA_df$tumourID, levels = summary_df$tumourID)
  toPlot$sample <- factor(toPlot$sample, levels = summary_df$tumourID)

  p1 <- ggplot(toPlot, aes(x = pipeline, y = sample, fill = count)) + 
    geom_tile() + 
    theme_minimal() +
    theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 17),
          legend.title = element_text(size = 10),
          legend.key.size = unit(3, "mm")) +
    scale_fill_viridis("log2(Counts)", limits=c(0, 13), option="mako", direction = -1) +
    #guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5, title = "log2(Counts + 1)")) +
    labs(x = "\nProtocol", y = "Tumour Sample", title = label)


  p2 <- ggplot() +
      geom_col(data = summary_df, aes(x = log2(mean_pct+1), y = tumourID), fill = "steelblue") +
      geom_point(data = rRNA_df, aes(x = log2(rRNA_Contamination_Pct+1), y = tumourID),
                  size = 1.5, color = "#4E598C", alpha = 0.8) +   
      labs(x = "log2+1 Normalized\nrRNA Contamination (%)", y = "Sample") +
      xlim(c(0, 2)) +
      theme_bw() + 
      theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
      )

  p <- p1 + p2 + plot_layout(width = c(1,2), guides = "collect")
  filename <- paste0("results/figures/suppfig5/rRNA_contamination_", label, ".png")
  ggsave(filename, p, width = 4, height = 5)
}


plot_rRNA_contamination(polyA_df, polyA_rRNA, "poly(A)-selection")
plot_rRNA_contamination(ribo0_df, ribo0_rRNA, "rRNA-depletion")

############################################################
# Plot seqdepth and counts
############################################################

toPlot <- polyA_df
lib_df <- polyA_lib
label = "poly(A)-selection"

plot_seqdepth <- function(toPlot, lib_df, label) {

  lib_df <- lib_df[order(lib_df$avg_counts),]
  lib_df$tumourID <- factor(lib_df$tumourID, levels = lib_df$tumourID)
  toPlot$sample <- factor(toPlot$sample, levels = lib_df$tumourID)

  p1 <- ggplot(toPlot, aes(x = pipeline, y = sample, fill = count)) + 
    geom_tile() + 
    theme_minimal() +
    theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 17),
          legend.title = element_text(size = 10),
          legend.key.size = unit(3, "mm")) +
    scale_fill_viridis("log2(Counts)", limits=c(0, 13), option="mako", direction = -1) +
    #guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5, title = "log2(Counts + 1)")) +
    labs(x = "\nProtocol", y = "Tumour Sample", title = label)


  p2 <- ggplot() +
      geom_col(data = lib_df, aes(x = avg_counts, y = tumourID), fill = "steelblue") +
      labs(x = "Average Read Counts", y = "Sample") +
      xlim(c(0, 78000000)) +
      theme_bw() + 
      theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8)
      )

  p <- p1 + p2 + plot_layout(width = c(1,2), guides = "collect")
  filename <- paste0("results/figures/suppfig5/seqdepth_", label, ".png")
  ggsave(filename, p, width = 4, height = 5)
}


plot_seqdepth(polyA_df, polyA_lib, "poly(A)-selection")
plot_seqdepth(ribo0_df, ribo0_lib, "rRNA-depletion")
