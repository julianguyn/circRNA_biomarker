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

valid <- c("cell_sub", "cell_all", "lungs")
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
  cell_sub = "../data/processed_cellline/common_samples",
  cell_all = "../data/processed_cellline/all_samples/",
  lungs = "../data/processed_lung"
)

files <- list.files(dir, recursive = TRUE, pattern = ".*counts.tsv")

############################################################
# Subset samples of interest if mapping all
############################################################

if (analysis == "cell_all") {
    # load in CIRI2 files
    ciri_gcsi <- fread(paste0(dir, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F) 
    ciri_gdsc <- fread(paste0(dir, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)
    ciri_ccle <- fread(paste0(dir, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)

    # find cell lines that overlap in at least 2 PSets
    gcsi_ccle <- intersect(ciri_gcsi$sample,ciri_ccle$sample)   # 481 samples
    gcsi_gdsc <- intersect(ciri_gcsi$sample,ciri_gdsc$sample)   # 112 samples
    ccle_gdsc <- intersect(ciri_ccle$sample,ciri_gdsc$sample)   # 91 samples

    to_keep <- unique(c(gcsi_ccle, gcsi_gdsc, ccle_gdsc))       # 588 total
    save(to_keep, file = "../data/temp/all_intersect_cells.RData")

    # remove files
    rm(ciri_gcsi, ciri_gdsc, ciri_ccle)
}

# helper function to filter for intersected cells
filter_circ <- function(counts, to_keep) {
    print("---Filtering for common cell lines")
    counts <- counts[counts$sample %in% to_keep, ]
    return(counts)
}

############################################################
# Process each file
############################################################

for (file in files) {

  print(paste("Starting file:", file))
  counts <- fread(paste0(dir, "/", file), data.table = FALSE)

  # filter for intersectef cells
  if (analysis == "cell_all") counts <- filter_circ(counts)

  # create GRanges of circRNA genomic coordinates
  gr <- IDstoBED(counts)

  # map circRNAs to transcript coordinates
  map <- circToGene(gr)

  # convert circRNA expression matrix labels to genes
  circExpToGene(ciri_gcsi_map, ciri_gcsi_sub, "CIRI2", "ciri_gcsi")
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

############################################################
# Convert circRNA expression matrix labels to genes
############################################################

# function to map circRNA expression matrices to gene-level annotations
circExpToGene <- function(circ_map, circ_exp) {

    print("---Making circRNA-mapped-gene expression matrix")

    # get dir labels
    dir <- sub("/.*", "", file) # should return pipeline or 'lung' TODO::check this and below
    label <- sub("_counts.tsv", "", sub(paste0(dir, "/"), "",  file)) # should return 'ciri_gcsi'

    # get circRNA-mapped-gene expression
    gene_exp <- circ_exp %>%
        pivot_longer(
            cols = -sample, 
            names_to = "old_name", 
            values_to = "expression"
        ) %>%
        inner_join(circ_map, by = "old_name", relationship = "many-to-many") %>%
        group_by(sample, gene) %>%
        summarize(expression = mean(expression), .groups = "drop") %>%
        pivot_wider(names_from = gene, values_from = expression) %>% 
        as.data.frame()

    # save output
    dir <- ifelse(dir != "lung", paste0("processed_cellline/GE_common_samples/", dir),"processed_lung/GE")
    filename <- paste0("../data/", dir, "/", label, "_counts.tsv")
    print(paste("Saving output to", filename))
    write.table(gene_exp, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}

# cells
circExpToGene(ciri_gcsi_map, ciri_gcsi_sub, "CIRI2", "ciri_gcsi")
circExpToGene(ciri_gdsc_map, ciri_gdsc_sub, "CIRI2", "ciri_gdsc")
circExpToGene(ciri_ccle_map, ciri_ccle_sub, "CIRI2", "ciri_ccle")

circExpToGene(circ_gcsi_map, circ_gcsi_sub, "CIRCexplorer2", "circ_gcsi")
circExpToGene(circ_gdsc_map, circ_gdsc_sub, "CIRCexplorer2", "circ_gdsc")
circExpToGene(circ_ccle_map, circ_ccle_sub, "CIRCexplorer2", "circ_ccle")

circExpToGene(cfnd_gcsi_map, cfnd_gcsi_sub, "circRNA_finder", "cfnd_gcsi")
circExpToGene(cfnd_gdsc_map, cfnd_gdsc_sub, "circRNA_finder", "cfnd_gdsc")
circExpToGene(cfnd_ccle_map, cfnd_ccle_sub, "circRNA_finder", "cfnd_ccle")

circExpToGene(fcrc_gcsi_map, fcrc_gcsi_sub, "find_circ", "fcrc_gcsi")
circExpToGene(fcrc_gdsc_map, fcrc_gdsc_sub, "find_circ", "fcrc_gdsc")
circExpToGene(fcrc_ccle_map, fcrc_ccle_sub, "find_circ", "fcrc_ccle")

# lung
circExpToGene(ciri_polyA_map, ciri_polyA, "lung", "ciri_polyA")
circExpToGene(ciri_ribo0_map, ciri_ribo0, "lung", "ciri_ribo0")

circExpToGene(circ_polyA_map, circ_polyA, "lung", "circ_polyA")
circExpToGene(circ_ribo0_map, circ_ribo0, "lung", "circ_ribo0")

circExpToGene(cfnd_polyA_map, cfnd_polyA, "lung", "cfnd_polyA")
circExpToGene(cfnd_ribo0_map, cfnd_ribo0, "lung", "cfnd_ribo0")

circExpToGene(fcrc_polyA_map, fcrc_polyA, "lung", "fcrc_polyA")
circExpToGene(fcrc_ribo0_map, fcrc_ribo0, "lung", "fcrc_ribo0")
