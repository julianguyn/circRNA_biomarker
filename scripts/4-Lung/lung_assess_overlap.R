# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(viridis)
    library(reshape2)
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

options(stringsAsFactors = FALSE)
#source("utils/palettes.R")

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

transcript_region_enrichment <- function(polyA_df, ribo0_df) {

  # get transcripts
  common <- intersect(colnames(polyA_df), colnames(ribo0_df))
  polyA_only <- colnames(polyA_df)[-which(colnames(polyA_df) %in% common)]
  ribo0_only <- colnames(ribo0_df)[-which(colnames(ribo0_df) %in% common)]

  # check genomic regions
  toPlot <- rbind(
    bin_enrichment(common, "Common transcripts"),
    bin_enrichment(polyA_only, "poly(A)-selection only"),
    bin_enrichment(ribo0_only, "rRNA-depleted only")
  )

  return(toPlot)
}

ciri <- transcript_region_enrichment(ciri_polyA, ciri_ribo0)
circ <- transcript_region_enrichment(circ_polyA, circ_ribo0)
cfnd <- transcript_region_enrichment(cfnd_polyA, cfnd_ribo0)
fcrc <- transcript_region_enrichment(fcrc_polyA, fcrc_ribo0)

save(ciri, circ, cfnd, fcrc, file = "assess_overlap.RData")