# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(reshape2)
    library(ggnewscale)
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(101)
source("utils/palettes.R")

analysis <- "circ"

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
# Helper function to assess bin enrichment
############################################################

# liftover from https://genome.ucsc.edu/cgi-bin/hgLiftOver
circbase <- read.table("../data/circBase/hglft_genome_1a58e3_fbb910.bed")
cb <- paste(circbase$V1, circbase$V2, circbase$V3, sep = ".")

bin_enrichment <- function(transcripts, label, bin_size = 1e6, n_perm = 10000) {

    parsed <- do.call(rbind, strsplit(transcripts, "\\."))

    gr <- GRanges(
        seqnames = parsed[, 1],
        ranges   = IRanges(start = as.numeric(parsed[, 2]),
                            end   = as.numeric(parsed[, 3]))
    )
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]

    all_circrnas <- GRanges(
      seqnames = circbase$V1,
      ranges = IRanges(start = circbase$V2, end = circbase$V3)
    )
    all_circrnas <- keepSeqlevels(all_circrnas, seqlevels(gr), pruning.mode = "coarse")

    bins <- tileGenome(seqlengths(gr), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)

    # count regions overlapping in bins
    obs_counts <- countOverlaps(bins, gr)
    n_regions <- length(gr)

    # permutation
    perm_counts <- matrix(NA, nrow = length(bins), ncol = n_perm)

    for (i in seq_len(n_perm)) {
      rand_idx <- sample(length(all_circrnas), n_regions, replace = FALSE)
      rand_gr  <- all_circrnas[rand_idx]
      perm_counts[, i] <- countOverlaps(bins, rand_gr)
      if (i %% 1000 == 0) message("Perm ", i, " of ", n_perm)
    }

    # empirical pval per bin 
    # p = (number of permutations with count >= observed, plus 1) / (n_perm + 1)
    pvals <- sapply(seq_along(bins), function(i) {
      (sum(perm_counts[i, ] >= obs_counts[i]) + 1) / (n_perm + 1)
    })
    padj <- p.adjust(pvals, method = "BH")

    # compile results
    results <- data.frame(
        chr = as.character(seqnames(bins)),
        start = start(bins),
        end = end(bins),
        count = obs_counts,
        pval = pvals,
        padj = padj,
        label = label
    )
    results$enriched <- ifelse(results$padj < 0.05 & results$count > 0, "Enriched", "Not Enriched")
    results <- results[order(results$padj),]
    return(results)
}

# v2: poisson distribution
bin_enrichment <- function(transcripts, label, bin_size = 1e6) {
    parsed <- do.call(rbind, strsplit(transcripts, "\\."))
    gr <- GRanges(
        seqnames = parsed[, 1],
        ranges   = IRanges(start = as.numeric(parsed[, 2]),
                            end   = as.numeric(parsed[, 3]))
    )

    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
    bins <- tileGenome(seqlengths(gr), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)

    # count regions overlapping in bins
    counts <- countOverlaps(bins, gr)

    # set up null model
    n_regions <- length(gr)
    total_bins <- length(bins)
    lambda <- n_regions / total_bins # expected count per bin (Poisson rate)

    # poisson test per bin
    pvals <- sapply(counts, function(x) {
        poisson.test(x, r = lambda, alternative = "greater")$p.value
    })
    padj <- p.adjust(pvals, method = "BH")

    # compile results
    results <- data.frame(
        chr = as.character(seqnames(bins)),
        start = start(bins),
        end = end(bins),
        count = counts,
        pval = pvals,
        padj = padj,
        label = label
    )
    results$enriched <- ifelse(results$padj < 0.05 & results$count > 0, "Enriched", "Not Enriched")
    results <- results[order(results$padj),]
    return(results)
}

############################################################
# Assess overlap between protocols
############################################################

transcript_region_enrichment_protocol <- function(polyA_df, ribo0_df) {

  # get transcripts
  common <- intersect(colnames(polyA_df), colnames(ribo0_df))
  polyA_only <- colnames(polyA_df)[-which(colnames(polyA_df) %in% common)]
  ribo0_only <- colnames(ribo0_df)[-which(colnames(ribo0_df) %in% common)]

  # check genomic regions
  toPlot <- rbind(
    bin_enrichment(common, "Both"),
    bin_enrichment(polyA_only, "poly(A)-selection only"),
    bin_enrichment(ribo0_only, "rRNA-depleted only")
  )

  return(toPlot)
}

message("%%%%%CIRI2")
ciri <- transcript_region_enrichment_protocol (ciri_polyA, ciri_ribo0)
message("%%%%%CIRCexplorer2")
circ <- transcript_region_enrichment_protocol (circ_polyA, circ_ribo0)
message("%%%%%circRNA_finder")
cfnd <- transcript_region_enrichment_protocol (cfnd_polyA, cfnd_ribo0)
message("%%%%%find_circ")
fcrc <- transcript_region_enrichment_protocol (fcrc_polyA, fcrc_ribo0)

save(ciri, circ, cfnd, fcrc, file = "../results/data/assess_overlap/overlap_protocol.RData")

############################################################
# Assess overlap between pipelines
############################################################

transcript_region_enrichment_pipeline <- function(ciri_df, circ_df, cfnd_df, fcrc_df) {

  # isolate transcripts
  ciri <- colnames(ciri_df)
  circ <- colnames(circ_df)
  cfnd <- colnames(cfnd_df)
  fcrc <- colnames(fcrc_df)

  # get common transcripts and remove from pool
  common <- intersect(intersect(ciri, circ), intersect(cfnd, fcrc))
  ciri <- ciri[-which(ciri %in% common)]
  circ <- circ[-which(circ %in% common)]
  cfnd <- cfnd[-which(cfnd %in% common)]
  fcrc <- fcrc[-which(fcrc %in% common)]

  # get transcripts in 3 pipelines
  ov1 <- intersect(ciri, intersect(circ, cfnd))
  ov2 <- intersect(ciri, intersect(circ, fcrc))
  ov3 <- intersect(ciri, intersect(cfnd, fcrc))
  ov4 <- intersect(circ, intersect(cfnd, fcrc))

  # check genomic regions
  toPlot <- rbind(
    bin_enrichment(common, "All Pipelines"),
    bin_enrichment(ov1, "Group1"),
    bin_enrichment(ov2, "Group2"),
    bin_enrichment(ov3, "Group3"),
    bin_enrichment(ov4, "Group4")
  )

  return(toPlot)
}

message("%%%%%polyA")
polyA <- transcript_region_enrichment_pipeline(ciri_polyA, circ_polyA, cfnd_polyA, fcrc_polyA)
message("%%%%%ribo0")
ribo0 <- transcript_region_enrichment_pipeline(ciri_ribo0, circ_ribo0, cfnd_ribo0, fcrc_ribo0)

save(polyA, ribo0, file = "../results/data/assess_overlap/overlap_pipeline.RData")

############################################################
# Plot Manhattan plots
############################################################

load("results/data/overlap_protocol.RData")
load("results/data/overlap_pipeline.RData")

plot_manhattan <- function(toPlot, pipeline, sig_threshold = 0.05) {

    toPlot$chr <- factor(toPlot$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
    toPlot <- toPlot[!is.na(toPlot$chr), ] 
    toPlot$padj <- ifelse(toPlot$padj < 1e-50, 1e-50, toPlot$padj)

    # get cumulative genomic positions for x axis
    data_cum <- toPlot %>%
      group_by(chr) %>%
      summarise(max_bp = as.numeric(max(end))) %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
      dplyr::select(chr, bp_add)

    toPlot <- toPlot %>%
      mutate(start = as.numeric(start)) %>%
      left_join(data_cum, by = "chr") %>%
      mutate(bp_cum = start + bp_add)

    # centre axis labels
    axis_set <- toPlot %>%
      group_by(chr) %>%
      summarise(center = mean(bp_cum))

    # try making plot
    p <- ggplot() +
      geom_point(data = toPlot, aes(x = bp_cum, y = -log10(padj), color = chr), alpha = 0.6) +
      scale_x_continuous(labels = axis_set$chr, breaks = axis_set$center) +
      scale_color_manual(values = rep(c("#CCCCCC", "#ABA8A8"), length(unique(toPlot$chr)))) +
      guides(color = "none") +
      geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "black") +
      new_scale_color() +
      geom_point(data = toPlot[toPlot$enriched == "Enriched",], aes(x = bp_cum, y = -log10(padj), color = label), alpha = 0.6) +
      scale_color_manual("", values = c(protocol_pal, pipeline_pal)) +
      labs(
        x = "Chromosome",
        y = expression(-log[10](FDR)),
        title = pipeline
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)
      )
    return(p)
}

# v2: plot counts
plot_manhattan <- function(toPlot, pipeline, sig_threshold = 0.05) {

    toPlot$chr <- factor(toPlot$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
    toPlot <- toPlot[!is.na(toPlot$chr), ] 

    # get cumulative genomic positions for x axis
    data_cum <- toPlot %>%
      group_by(chr) %>%
      summarise(max_bp = as.numeric(max(end))) %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
      dplyr::select(chr, bp_add)

    toPlot <- toPlot %>%
      mutate(start = as.numeric(start)) %>%
      left_join(data_cum, by = "chr") %>%
      mutate(bp_cum = start + bp_add)

    # centre axis labels
    axis_set <- toPlot %>%
      group_by(chr) %>%
      summarise(center = mean(bp_cum))

    # try making plot
    p <- ggplot() +
      geom_point(data = toPlot, aes(x = bp_cum, y = log2(count), color = chr), alpha = 0.6) +
      scale_x_continuous(labels = axis_set$chr, breaks = axis_set$center) +
      scale_color_manual(values = rep(c("#CCCCCC", "#ABA8A8"), length(unique(toPlot$chr)))) +
      guides(color = "none") +
      geom_hline(yintercept = log2(5), linetype = "dashed", color = "black") +
      new_scale_color() +
      geom_point(data = toPlot[which(log2(toPlot$count) > log2(5) & toPlot$label == "rRNA-depleted only"),], aes(x = bp_cum, y = log2(count), color = label), alpha = 0.6) +
      geom_point(data = toPlot[which(log2(toPlot$count) > log2(5) & toPlot$label != "rRNA-depleted only"),], aes(x = bp_cum, y = log2(count), color = label), alpha = 0.6) +
      scale_color_manual("", values = c(protocol_pal, pipeline_pal)) +
      labs(
        x = "Chromosome",
        y = expression(-log[2](Count)),
        title = pipeline
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)
      )
    return(p)
}

p1 <- plot_manhattan(ciri, "CIRI2")
p2 <- plot_manhattan(circ, "CIRCexplorer2")
p3 <- plot_manhattan(cfnd, "circRNA_finder")
p4 <- plot_manhattan(fcrc, "find_circ")

p5 <- plot_manhattan(polyA, "polyA")
p6 <- plot_manhattan(ribo0, "ribo0")

p1 <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")
p2 <- (p5 / p6)+ plot_layout(guides = "collect")

filename <- "results/figures/suppfig5/overlap_protocol-count.png"
ggsave(filename, p1, w=12, h=4)

filename <- "results/figures/suppfig5/overlap_pipeline-count.png"
ggsave(filename, p2, w=6.5, h=4)
