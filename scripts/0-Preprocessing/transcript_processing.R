# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
})

############################################################
# Load data
############################################################

# load in common 48 cell lines
load("../data/processed_cellline/common_samples/inter.RData")

# read in PSets
gcsi <- readRDS("../data/PSets/gCSI.rds") |> updateObject()
ccle <- readRDS("../data/PSets/CCLE.rds") |> updateObject()
gdsc <- readRDS("../data/PSets/GDSC2-8.2.rds") |> updateObject()

############################################################
# Get molecular profiles
############################################################

# summarize isoform expression
gcsi_rna <- summarizeMolecularProfiles(gcsi, mDataType = "Kallisto_0.46.1.isoforms.counts", cell.lines = intersected_rnacells)
gdsc_rna <- summarizeMolecularProfiles(ccle, mDataType = "Kallisto_0.46.1.isoforms.counts", cell.lines = intersected_rnacells)
gdsc_rna <- summarizeMolecularProfiles(gdsc, mDataType = "Kallisto_0.46.1.isoforms.counts", cell.lines = intersected_rnacells)

# summarize gene expression
gcsi_gene <- summarizeMolecularProfiles(gcsi, mDataType = "Kallisto_0.46.1.rnaseq.counts", cell.lines = intersected_rnacells)
ccle_gene <- summarizeMolecularProfiles(ccle, mDataType = "Kallisto_0.46.1.rnaseq.counts", cell.lines = intersected_rnacells)
gdsc_gene <- summarizeMolecularProfiles(gdsc, mDataType = "Kallisto_0.46.1.rnaseq.counts", cell.lines = intersected_rnacells)


############################################################
# Keep only protein-coding transcripts
############################################################

# save corresponding gene_id for gene expression quantification
transcript_gene <- gdsc_rna@elementMetadata[,c("gene_id", "transcript_id")]
protein_transcripts <- gcsi_rna@elementMetadata$transcript_id[gCSI_rna@elementMetadata$transcript_type == "protein_coding"]
protein_genes <- gdsc_rna@elementMetadata$gene_id[which(gdsc_rna@elementMetadata$gene_type == "protein_coding" & gdsc_rna@elementMetadata$gene_id %in% transcript_gene$gene_id)]

# get isoform transcripts
expr_gcsi_i <- gcsi_rna@assays@data$expr[protein_transcripts,]
expr_ccle_i <- ccle_rna@assays@data$expr[protein_transcripts,]
expr_gdsc_i <- gdsc_rna@assays@data$expr[protein_transcripts,]

# get gene transcripts
expr_gcsi_p <- gcsi_gene@assays@data$expr[protein_genes,]
expr_ccle_p <- ccle_gene@assays@data$expr[protein_genes,]
expr_gdsc_p <- gdsc_gene@assays@data$expr[protein_genes,]

############################################################
# Filter for low expression
############################################################

# remove transcripts with 0 expression across all cell lines
expr_gcsi_i <- expr_gcsi_i[apply(expr_gcsi_i[,-1], 1, function(x) !all(x==0)),]
expr_ccle_i <- expr_ccle_i[apply(expr_ccle_i[,-1], 1, function(x) !all(x==0)),]
expr_gdsc_i <- expr_gdsc_i[apply(expr_gdsc_i[,-1], 1, function(x) !all(x==0)),]

expr_gcsi_p <- expr_gcsi_p[apply(expr_gcsi_p[,-1], 1, function(x) !all(x==0)),]
expr_ccle_p <- expr_ccle_p[apply(expr_ccle_p[,-1], 1, function(x) !all(x==0)),]
expr_gdsc_p <- expr_gdsc_p[apply(expr_gdsc_p[,-1], 1, function(x) !all(x==0)),]

# only keep transcripts found in all datasets 
intersected_t <- Reduce(intersect, list(rownames(expr_gcsi_i), rownames(expr_ccle_i), rownames(expr_gdsc_i)))
intersected_g <- Reduce(intersect, list(rownames(expr_gcsi_p), rownames(expr_ccle_p), rownames(expr_gdsc_p)))

expr_gcsi_i <- expr_gcsi_i[intersected_t,]
expr_ccle_i <- expr_ccle_i[intersected_t,]
expr_gdsc_i <- expr_gdsc_i[intersected_t,]

expr_gcsi_p <- expr_gcsi_p[intersected_g,]
expr_ccle_p <- expr_ccle_p[intersected_g,]
expr_gdsc_p <- expr_gdsc_p[intersected_g,]

############################################################
# Remove housekeeping transcripts from isoforms
############################################################

# remove housekeeping transcripts (defined by HRT Atlas)
transcripts_noid <- gsub("\\..*","",intersected_t)
housekeeping <- read.csv("../data/transcript_stability/Housekeeping_TranscriptsHuman.csv", sep = ";")

expr_gcsi_i <- expr_gcsi_i[which(!transcripts_noid %in%  housekeeping$Ensembl),]
expr_ccle_i <- expr_ccle_i[which(!transcripts_noid %in%  housekeeping$Ensembl),]
expr_gdsc_i <- expr_gdsc_i[which(!transcripts_noid %in%  housekeeping$Ensembl),]


############################################################
# Format for stability index calculations
############################################################

expr_gcsi_i <- t(expr_gcsi_i) |> as.data.frame()
expr_ccle_i <- t(expr_ccle_i) |> as.data.frame()
expr_gdsc_i <- t(expr_gdsc_i) |> as.data.frame()

expr_gcsi_p <- t(expr_gcsi_p) |> as.data.frame()
expr_ccle_p <- t(expr_ccle_p) |> as.data.frame()
expr_gdsc_p <- t(expr_gdsc_p) |> as.data.frame()

############################################################
# Save matrices
############################################################

save(expr_gcsi_i, expr_ccle_i, expr_gdsc_i, file = "../results/data/isoform_expression.RData")
save(expr_gcsi_p, expr_ccle_p, expr_gdsc_p, file = "../results/data/gene_expression.RData")


############################################################
# Save metadata (for feature influence)
############################################################

isoform_meta <- gcsi_rna@elementMetadata
gene_meta <- gcsi_gene@elementMetadata

save(isoform_meta, gene_meta, file = "../data/rnaseq_meta/transcript_meta.RData")