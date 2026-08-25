# evaluate influence of sequencing depth on results

suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
})

options(stringsAsFactors = FALSE)

source("utils/circ_processing.R")

analysis <- "cells"
dir <- "../data/raw_cellline"

############################################################
# Load in cell line metadata
############################################################

load("../results/data/tissue-metadata.RData")
load("../data/processed_cellline/common_samples/inter.RData")

############################################################
# Load in metadata
############################################################

# get input dirs
files <- list.files(dir, recursive = TRUE, pattern = ".*counts.tsv")

for (file in files) {

  print(paste("Starting file:", file))

  # extract variables
  pipeline <- sub("/.*", "", file)
  count_file <- sub(".*/", "", file)
  dataset <- sub(".*_", "", sub(".*/", "", sub("_counts.tsv", "", file)))

  counts <- fread(paste0(dir, "/", file), data.table = FALSE)
  colnames(counts)[1] <- "sample"

  # get library sizes for dataset
  lib <- switch(
    dataset,
    gcsi = "../data/raw_cellline/readCounts/gcsi.tsv",
    ccle = "../data/raw_cellline/readCounts/ccle.tsv",
    gdsc = "../data/raw_cellline/readCounts/gdsc.tsv"
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
    reads$sample <- gdsc$cellid[match(reads$sample, sub("#.*", "", rownames(gdsc)))]
  }

  # subset for intersected cell lines
  counts <- counts[counts$sample %in% intersected_rnacells,]
  reads <- reads[reads$sample %in% intersected_rnacells,]

  # average across technical replicates
  counts <- avg_reps(counts)

  # save raw counts counts
  outdir <- sub("raw", "seqdepth", dir)
  if (!dir.exists(paste0(outdir, "/", pipeline))) {
    dir.create(paste0(outdir, "/", pipeline), recursive = TRUE)
  }
  outfile <- paste0(outdir, "/", pipeline, "/", count_file)
  print(paste("Saving raw counts to:", outfile))
  write.table(
    counts,
    file = outfile,
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
  readfile <- paste0(outdir, sub(".*raw_cellline", "", lib))
  write.table(
    reads,
    file = readfile,
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
}
print("done")