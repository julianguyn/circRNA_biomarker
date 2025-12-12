# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(stringr)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
    library(tidyr)
})

source("utils/circ_to_gene.R")

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("cells", "lungs")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

############################################################
# Load in genomic annotation data
############################################################

# load in hg38 genomic coordinates
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")

# map gene symbols
mcols(exons)$gene_name <- mapIds(
    org.Hs.eg.db,
    keys = names(exons),
    column = "SYMBOL",
    keytype = "ENTREZID"
)

############################################################
# Locate input files
############################################################

# get input dirs
dir <- switch(
  analysis,
  cells = "../data/processed_cellline/common_samples",
  lungs = "../data/processed_lung"
)

files <- list.files(dir, recursive = TRUE, pattern = ".*counts.tsv")

############################################################
# Process each file
############################################################

for (file in files) {

  print(paste("Starting file:", file))
  counts <- fread(paste0(dir, "/", file), data.table = FALSE)

  # create GRanges of circRNA genomic coordinates
  gr <- IDstoBED(counts)

  # map circRNAs to transcript coordinates
  map <- circToGene(gr)

  # convert circRNA expression matrix labels to genes
  circExpToGene(map, counts)
}

############################################################
# Saved map checkpoints (from time out)
############################################################

# cell_all
#save(ciri_gcsi_map, ciri_gdsc_map, ciri_ccle_map,
#     circ_gcsi_map, circ_gdsc_map, circ_ccle_map,
#     cfnd_gcsi_map, cfnd_gdsc_map, cfnd_ccle_map,
#     fcrc_gcsi_map, fcrc_gdsc_map, fcrc_ccle_map,
#     file = "../results/data/temp/circ_to_gene_map.RData")

# lungs
#save(ciri_polyA_map, ciri_ribo0_map,
#     circ_polyA_map, circ_ribo0_map,
#     cfnd_polyA_map, cfnd_ribo0_map,
#     fcrc_polyA_map, fcrc_ribo0_map,
#     file = "../results/data/temp/circ_lung_to_gene_map.RData")
