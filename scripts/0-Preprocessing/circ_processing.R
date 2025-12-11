#' Create analysis-ready count matrices for circRNA
#' 
suppressPackageStartupMessages({
    library(data.table)
})

options(stringsAsFactors = FALSE)

source("utils/circ_processing.R")

# load libraries
#suppressMessages(library(data.table))
#suppressMessages(library(ggplot2))s
#suppressMessages(library(ggpubr))
#suppressMessages(library(PharmacoGx))
#suppressMessages(library(stringr))
#suppressMessages(library(ComplexHeatmap))

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
# Load in metadata
############################################################

#gcsi <- readRDS("gCSI.rds") |> updateObject()
#gcsi_metadata <- gcsi@molecularProfiles$Kallisto_0.46.1.rnaseq@colData 
#write.table(gcsi_metadata, file = "../data/rnaseq_meta/tissue_meta/gcsi_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#gdsc <- readRDS("GDSC2-8.2.rds") |> updateObject()
#gdsc_metadata <- gdsc@molecularProfiles$Kallisto_0.46.1.rnaseq@colData
#write.table(gdsc_metadata, file = "../data/rnaseq_meta/tissue_meta/gdsc_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#ccle <- readRDS("CCLE.rds") |> updateObject()
#ccle_metadata <- ccle@molecularProfiles$Kallisto_0.46.1.rnaseq@colData
#write.table(ccle_metadata, file = "../data/rnaseq_meta/tissue_meta/ccle_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)


# load in metadata for tissue-specific analysis
gcsi_tmeta <- fread("../data/rnaseq_meta/tissue_meta/gcsi_metadata.tsv")
ccle_tmeta <- fread("../data/rnaseq_meta/tissue_meta/ccle_metadata.tsv")
gdsc_tmeta <- fread("../data/rnaseq_meta/tissue_meta/gdsc_metadata.tsv")


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
  counts <- fread(paste0(dir, "/", file), data.table = FALSE)
  colnames(counts)[1] <- "sample"

  # catch the ribo0 lungs
  counts$sample <- sub("_INPUT", "", counts$sample)

  # extract variables
  pipeline <- sub("/.*", "", file)
  count_file <- sub(".*/", "", file)
  dataset <- sub(".*_", "", sub(".*/", "", sub("_counts.tsv", "", file)))
  label <- paste(analysis, dataset, sep = "_")

  # get library sizes for dataset
  lib <- switch(
    label,
    cells_gcsi = "../data/raw_cellline/readCounts/gcsi.tsv",
    cells_ccle = "../data/raw_cellline/readCounts/ccle.tsv",
    cells_gdsc = "../data/raw_cellline/readCounts/gdsc.tsv",
    lungs_polyA = "../data/raw_lung/readCounts/polyA.tsv",
    lungs_ribo0 = "../data/raw_lung/readCounts/ribo0.tsv"
  )
  print(paste("Reading in library file:", lib))
  reads <- read.table(lib, header = TRUE)

  # normalize by read depth
  counts <- get_cpm(counts, reads)

  # average across technical replicates
  counts <- avg_reps(counts)

  # save normalized counts
  outdir <- sub("raw", "processed", dir)
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
}
print("done")



############################################################
# Load in data
############################################################

# read in circRNA matrices
suppressWarnings(ciri_gcsi <- fread("../data/raw_cellline/CIRI2/ciri_gcsi_counts.tsv", data.table = F))
suppressWarnings(ciri_gdsc <- fread("../data/raw_cellline/CIRI2/ciri_gdsc_counts.tsv", data.table = F))
suppressWarnings(ciri_ccle <- fread("../data/raw_cellline/CIRI2/ciri_ccle_counts.tsv", data.table = F))

suppressWarnings(circ_gcsi <- fread("../data/raw_cellline/CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F))
suppressWarnings(circ_gdsc <- fread("../data/raw_cellline/CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F))
suppressWarnings(circ_ccle <- fread("../data/raw_cellline/CIRCexplorer2/circ_ccle_counts.tsv", data.table = F))

suppressWarnings(cfnd_gcsi <- fread("../data/raw_cellline/circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F))
suppressWarnings(cfnd_gdsc <- fread("../data/raw_cellline/circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F))
suppressWarnings(cfnd_ccle <- fread("../data/raw_cellline/circRNA_finder/cfnd_ccle_counts.tsv", data.table = F))

suppressWarnings(fcrc_gcsi <- fread("../data/raw_cellline/find_circ/fcrc_gcsi_counts.tsv", data.table = F))
suppressWarnings(fcrc_gdsc <- fread("../data/raw_cellline/find_circ/fcrc_gdsc_counts.tsv", data.table = F))
suppressWarnings(fcrc_ccle <- fread("../data/raw_cellline/find_circ/fcrc_ccle_counts.tsv", data.table = F))


# load in metadata for tissue-specific analysis



# function to match cell.id to unique.cellid from PharmacoGx
matchToIDTable <- function(
  ids,
  tbl,
  column,
  returnColumn = "unique.cellid"
) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}

# ========== Normalize by read depth ========== #

# load in library counts
gcsi_reads <- read.table("/cluster/home/julian/circRNA_biomarker/data/raw_cellline/readCounts/gcsi.tsv")
ccle_reads <- read.table("/cluster/home/julian/circRNA_biomarker/data/raw_cellline/readCounts/ccle.tsv")
gdsc_reads <- read.table("/cluster/home/julian/circRNA_biomarker/data/raw_cellline/readCounts/gdsc.tsv")

# match sample names
gcsi_reads$sample <- gcsi$cellid[match(gcsi_reads$sample, rownames(gcsi))]
ccle_reads$sample <- ccle$cellid[match(ccle_reads$sample, rownames(ccle))]
gdsc_reads$sample <- gdsc$cellid[match(gdsc_reads$sample, rownames(gdsc))]



# normalize circ data frames
ciri_gcsi <- cpm(ciri_gcsi, gcsi_reads)
ciri_ccle <- cpm(ciri_ccle, ccle_reads)
ciri_gdsc <- cpm(ciri_gdsc, gdsc_reads)

circ_gcsi <- cpm(circ_gcsi, gcsi_reads)
circ_ccle <- cpm(circ_ccle, ccle_reads)
circ_gdsc <- cpm(circ_gdsc, gdsc_reads)

cfnd_gcsi <- cpm(cfnd_gcsi, gcsi_reads)
cfnd_ccle <- cpm(cfnd_ccle, ccle_reads)
cfnd_gdsc <- cpm(cfnd_gdsc, gdsc_reads)

fcrc_gcsi <- cpm(fcrc_gcsi, gcsi_reads)
fcrc_ccle <- cpm(fcrc_ccle, ccle_reads)
fcrc_gdsc <- cpm(fcrc_gdsc, gdsc_reads)


# ========== Average across all technical replicates ========== #

# average technical replicates across circ data frames
ciri_gcsi <- avg_reps(ciri_gcsi)
ciri_ccle <- avg_reps(ciri_ccle)
ciri_gdsc <- avg_reps(ciri_gdsc)

circ_gcsi <- avg_reps(circ_gcsi)
circ_ccle <- avg_reps(circ_ccle)
circ_gdsc <- avg_reps(circ_gdsc)

cfnd_gcsi <- avg_reps(cfnd_gcsi)
cfnd_ccle <- avg_reps(cfnd_ccle)
cfnd_gdsc <- avg_reps(cfnd_gdsc)

fcrc_gcsi <- avg_reps(fcrc_gcsi)
fcrc_ccle <- avg_reps(fcrc_ccle)
fcrc_gdsc <- avg_reps(fcrc_gdsc)


# save dataframes
write.table(ciri_gcsi, file = "../data/processed_cellline/all_samples/CIRI2/ciri_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_gdsc, file = "../data/processed_cellline/all_samples/CIRI2/ciri_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_ccle, file = "../data/processed_cellline/all_samples/CIRI2/ciri_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(circ_gcsi, file = "../data/processed_cellline/all_samples/CIRCexplorer2/circ_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_gdsc, file = "../data/processed_cellline/all_samples/CIRCexplorer2/circ_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_ccle, file = "../data/processed_cellline/all_samples/CIRCexplorer2/circ_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(cfnd_gcsi, file = "../data/processed_cellline/all_samples/circRNA_finder/cfnd_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_gdsc, file = "../data/processed_cellline/all_samples/circRNA_finder/cfnd_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_ccle, file = "../data/processed_cellline/all_samples/circRNA_finder/cfnd_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(fcrc_gcsi, file = "../data/processed_cellline/all_samples/find_circ/fcrc_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_gdsc, file = "../data/processed_cellline/all_samples/find_circ/fcrc_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_ccle, file = "../data/processed_cellline/all_samples/find_circ/fcrc_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)



# ========== Create dataframes of intersected cell lines ========== #

# get common cell lines across gCSI, CCLE, and GDSC
intersected_rnacells <- intersect(intersect(ciri_gcsi$sample, ciri_ccle$sample), ciri_gdsc$sample)
save(intersected_rnacells, file="../data/processed_cellline/common_samples/inter.RData")


# subset dataframes to keep only common cell lines
ciri_gcsi_sub <- ciri_gcsi[ciri_gcsi$sample %in% intersected_rnacells,]
ciri_ccle_sub <- ciri_ccle[ciri_ccle$sample %in% intersected_rnacells,]
ciri_gdsc_sub <- ciri_gdsc[ciri_gdsc$sample %in% intersected_rnacells,]

circ_gcsi_sub <- circ_gcsi[circ_gcsi$sample %in% intersected_rnacells,]
circ_ccle_sub <- circ_ccle[circ_ccle$sample %in% intersected_rnacells,]
circ_gdsc_sub <- circ_gdsc[circ_gdsc$sample %in% intersected_rnacells,]

cfnd_gcsi_sub <- cfnd_gcsi[cfnd_gcsi$sample %in% intersected_rnacells,]
cfnd_ccle_sub <- cfnd_ccle[cfnd_ccle$sample %in% intersected_rnacells,]
cfnd_gdsc_sub <- cfnd_gdsc[cfnd_gdsc$sample %in% intersected_rnacells,]

fcrc_gcsi_sub <- fcrc_gcsi[fcrc_gcsi$sample %in% intersected_rnacells,]
fcrc_ccle_sub <- fcrc_ccle[fcrc_ccle$sample %in% intersected_rnacells,]
fcrc_gdsc_sub <- fcrc_gdsc[fcrc_gdsc$sample %in% intersected_rnacells,]


### TODO: remove
ciri_gdsc_sub[is.na(ciri_gdsc_sub)] <- 0
circ_gdsc_sub[is.na(circ_gdsc_sub)] <- 0
cfnd_gdsc_sub[is.na(cfnd_gdsc_sub)] <- 0
fcrc_gdsc_sub[is.na(fcrc_gdsc_sub)] <- 0


# function to filter transcripts with 0 expression across intersected cell lines
filter_transcripts <- function(df) {

  # save cell line labels
  sample <- df$sample

  # remove cell line labels
  df$sample <- NULL

  # filter transcripts with 0 expression
  df_sub <- df[,-which(colnames(df) %in% names(which(colSums(df) == 0)))]

  # re-add cell line labels
  df <- cbind(sample, df_sub)

  return(df)
}

ciri_gcsi_sub <- filter_transcripts(ciri_gcsi_sub)
ciri_ccle_sub <- filter_transcripts(ciri_ccle_sub)
ciri_gdsc_sub <- filter_transcripts(ciri_gdsc_sub)

circ_gcsi_sub <- filter_transcripts(circ_gcsi_sub)
circ_ccle_sub <- filter_transcripts(circ_ccle_sub)
circ_gdsc_sub <- filter_transcripts(circ_gdsc_sub)

cfnd_gcsi_sub <- filter_transcripts(cfnd_gcsi_sub)
cfnd_ccle_sub <- filter_transcripts(cfnd_ccle_sub)
cfnd_gdsc_sub <- filter_transcripts(cfnd_gdsc_sub)

fcrc_gcsi_sub <- filter_transcripts(fcrc_gcsi_sub)
fcrc_ccle_sub <- filter_transcripts(fcrc_ccle_sub)
fcrc_gdsc_sub <- filter_transcripts(fcrc_gdsc_sub)


# save dataframes
write.table(ciri_gcsi_sub, file = "../data/processed_cellline/common_samples/CIRI2/ciri_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_gdsc_sub, file = "../data/processed_cellline/common_samples/CIRI2/ciri_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_ccle_sub, file = "../data/processed_cellline/common_samples/CIRI2/ciri_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(circ_gcsi_sub, file = "../data/processed_cellline/common_samples/CIRCexplorer2/circ_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_gdsc_sub, file = "../data/processed_cellline/common_samples/CIRCexplorer2/circ_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_ccle_sub, file = "../data/processed_cellline/common_samples/CIRCexplorer2/circ_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(cfnd_gcsi_sub, file = "../data/processed_cellline/common_samples/circRNA_finder/cfnd_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_gdsc_sub, file = "../data/processed_cellline/common_samples/circRNA_finder/cfnd_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_ccle_sub, file = "../data/processed_cellline/common_samples/circRNA_finder/cfnd_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(fcrc_gcsi_sub, file = "../data/processed_cellline/common_samples/find_circ/fcrc_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_gdsc_sub, file = "../data/processed_cellline/common_samples/find_circ/fcrc_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_ccle_sub, file = "../data/processed_cellline/common_samples/find_circ/fcrc_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
