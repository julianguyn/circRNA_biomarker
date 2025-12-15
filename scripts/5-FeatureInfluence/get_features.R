# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(stringr)
    library(BSgenome)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
    library(tidyr)
    library(PharmacoGx)
})

set.seed(200)

source("utils/get_features.R")

############################################################
# Load in data
############################################################

# load circRNA expression data (subsetted in si_distribution.R)
load("../results/data/temp/circ_stability_subsetdf2.RData")

# load gene and isoform expression
load("../results/data/isoform_expression.RData")
load("../results/data/gene_expression.RData")

# load in GRCh38 reference genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# load in hg38 genomic coordinates
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene") 

# map gene symbols
mcols(exons)$gene_name <- mapIds(org.Hs.eg.db, 
                                    keys = names(exons), 
                                    column = "SYMBOL", 
                                    keytype = "ENTREZID")

############################################################
# Get gene and isoform coords
############################################################

# load in gene and isoform metadata (from transcript_processing.R)
load("../data/rnaseq_meta/transcript_meta.RData")

# get coordinates
gene_meta$coords <- paste0(gene_meta$seqnames, ".", gene_meta$start, ".", gene_meta$end)
isoform_meta$coords <- paste0(isoform_meta$seqnames, ".", isoform_meta$start, ".", isoform_meta$end)

# set coords as colnames
colnames(expr_gcsi_p) <- gene_meta$coords[match(colnames(expr_gcsi_p), gene_meta$gene_id)]
colnames(expr_ccle_p) <- gene_meta$coords[match(colnames(ccle_gcsi_p), gene_meta$gene_id)]
colnames(expr_gdsc_p) <- gene_meta$coords[match(colnames(expr_gdsc_p), gene_meta$gene_id)]

colnames(expr_gcsi_i) <- gene_meta$coords[match(colnames(expr_gcsi_i), gene_meta$transcript_id_id)]
colnames(expr_ccle_i) <- gene_meta$coords[match(colnames(ccle_gcsi_i), gene_meta$transcript_id_id)]
colnames(expr_gdsc_i) <- gene_meta$coords[match(colnames(expr_gdsc_i), gene_meta$transcript_id_id)]

############################################################
# Get transcript and sequencing features
############################################################

# circRNA
ciri_gcsi_ft <- get_features(ciri_gcsi)
ciri_gdsc_ft <- get_features(ciri_gdsc)
ciri_ccle_ft <- get_features(ciri_ccle)

circ_gcsi_ft <- get_features(circ_gcsi)
circ_gdsc_ft <- get_features(circ_gdsc)
circ_ccle_ft <- get_features(circ_ccle)

cfnd_gcsi_ft <- get_features(cfnd_gcsi)
cfnd_gdsc_ft <- get_features(cfnd_gdsc)
cfnd_ccle_ft <- get_features(cfnd_ccle)

fcrc_gcsi_ft <- get_features(fcrc_gcsi)
fcrc_gdsc_ft <- get_features(fcrc_gdsc)
fcrc_ccle_ft <- get_features(fcrc_ccle)

# gene expression
gene_gcsi_ft <- get_features(expr_gcsi_p)
gene_gdsc_ft <- get_features(expr_gdsc_p)
gene_ccle_ft <- get_features(expr_ccle_p)

# isoform expression
isof_gcsi_ft <- get_features(expr_gcsi_i)
isof_gdsc_ft <- get_features(expr_gdsc_i)
isof_ccle_ft <- get_features(expr_ccle_i)

############################################################
# Create dataframes for Elasticnet
############################################################

# load the stability indices
load("../results/data/stability/gene_isoform_stability.RData")
load("../results/data/stability/circ_stability.RData")

format_stability(ciri_stability, ciri_gcsi_ft, ciri_ccle_ft, ciri_gdsc_ft, "ciri")
format_stability(circ_stability, circ_gcsi_ft, circ_ccle_ft, circ_gdsc_ft, "circ")
format_stability(cfnd_stability, cfnd_gcsi_ft, cfnd_ccle_ft, cfnd_gdsc_ft, "cfnd")
format_stability(fcrc_stability, fcrc_gcsi_ft, fcrc_ccle_ft, fcrc_gdsc_ft, "fcrc")

format_stability(gene_stability, gene_gcsi_ft, gene_ccle_ft, gene_gdsc_ft, "gene")
format_stability(transcript_stability, isof_gcsi_ft, isof_ccle_ft, isof_gdsc_ft, "transcript")