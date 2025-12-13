#' Create analysis-ready count matrices for circRNA
#' 
suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
    #library(PharmacoGx)
})

options(stringsAsFactors = FALSE)

source("utils/circ_processing.R")

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("cells", "lung")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}


############################################################
# Load in cell line metadata
############################################################

#gcsi <- readRDS("../data/PSets/gCSI.rds") |> updateObject()
#gcsi_metadata <- gcsi@molecularProfiles$Kallisto_0.46.1.rnaseq@colData 
#write.table(gcsi_metadata, file = "../data/rnaseq_meta/tissue_meta/gcsi_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#gdsc <- readRDS("../data/PSets/GDSC2-8.2.rds") |> updateObject()
#gdsc_metadata <- gdsc@molecularProfiles$Kallisto_0.46.1.rnaseq@colData
#write.table(gdsc_metadata, file = "../data/rnaseq_meta/tissue_meta/gdsc_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#ccle <- readRDS("../data/PSets/CCLE.rds") |> updateObject()
#ccle_metadata <- ccle@molecularProfiles$Kallisto_0.46.1.rnaseq@colData
#write.table(ccle_metadata, file = "../data/rnaseq_meta/tissue_meta/ccle_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

if (analysis == "cells") {
  # read in cell annotations
  gcsi <- load_anno("../data/rnaseq_meta/gcsi_rnaseq_meta.csv", "gcsi")
  ccle <- load_anno("../data/rnaseq_meta/ccle_rnaseq_meta.csv", "ccle")
  gdsc <- load_anno("../data/rnaseq_meta/gdsc/sample_file.csv", "gdsc")

  # save tissue metadata
  save(gcsi, ccle, gdsc, file = "../results/data/tissue-metadata.RData")
}

############################################################
# Load in intersected cell lines (n=48)
############################################################

# intersected_rnacells <- intersect(intersect(ciri_gcsi$sample, ciri_ccle$sample), ciri_gdsc$sample)
# save(intersected_rnacells, file="../data/processed_cellline/common_samples/inter.RData")
load("../data/processed_cellline/common_samples/inter.RData")

############################################################
# Load in metadata
############################################################

# get input dirs
dir <- switch(
  analysis,
  cells = "../data/raw_cellline",
  lungs = "../data/raw_lung"
)

files <- list.files(dir, recursive = TRUE, pattern = ".*counts.tsv")

for (file in files) {

  print(paste("Starting file:", file))

  # extract variables
  pipeline <- sub("/.*", "", file)
  count_file <- sub(".*/", "", file)
  dataset <- sub(".*_", "", sub(".*/", "", sub("_counts.tsv", "", file)))

  counts <- fread(paste0(dir, "/", file), data.table = FALSE)
  colnames(counts)[1] <- "sample"

  # catch the ribo0 lungs
  if (analysis == lungs) counts$sample <- sub("_INPUT", "", counts$sample)

  # get library sizes for dataset
  lib <- switch(
    dataset,
    gcsi = "../data/raw_cellline/readCounts/gcsi.tsv",
    ccle = "../data/raw_cellline/readCounts/ccle.tsv",
    gdsc = "../data/raw_cellline/readCounts/gdsc.tsv",
    polyA = "../data/raw_lung/readCounts/polyA.tsv",
    ribo0 = "../data/raw_lung/readCounts/ribo0.tsv"
  )
  print(paste("Reading in library file:", lib))
  reads <- read.table(lib, header = TRUE)

  # process cells
  if (dataset == "gcsi") {
    counts$sample <- gcsi$cellid[match(gsub("gcsi", "", counts$sample), rownames(gcsi))]
    reads$sample <- gcsi$cellid[match(reads$sample, rownames(gcsi))]
  } else if (dataset == "ccle") {
    counts$sample <- ccle$cellid[match(counts$sample, rownames(ccle))]
    reads$sample <- ccle$cellid[match(reads$sample, rownames(ccle))]
  } else if (dataset == "gdsc") {
    counts$sample <- gdsc$cellid[match(counts$sample, rownames(gdsc))]
    reads$sample <- gdsc$cellid[match(reads$sample, rownames(gdsc))]
  }

  # normalize by read depth
  counts <- get_cpm(counts, reads)

  # average across technical replicates
  counts <- avg_reps(counts)

  # save normalized counts
  outdir <- sub("raw", "processed", dir)
  if (analysis == "cells") outdir <- paste0(outdir, "/all_samples")
  if (!dir.exists(paste0(outdir, "/", pipeline))) {
    dir.create(paste0(outdir, "/", pipeline), recursive = TRUE)
  }
  outfile <- paste0(outdir, "/", pipeline, "/", count_file)
  print(paste("Saving normalized counts to:", outfile))
  write.table(
    counts,
    file = outfile,
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )

  # get intersected cell lines
  if (analysis == "cells") {
    counts <- counts[counts$sample %in% intersected_rnacells, ]
    counts[is.na(counts)] <- 0
    counts <- filter_transcripts(counts)
    outfile <- sub("all", "common", outfile)
    print(paste("Saving subsetted normalized counts to:", outfile))
    write.table(
      counts,
      file = outfile,
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
  }
  
}
print("done")