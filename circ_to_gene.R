

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


############################################################
# Load in data 
############################################################

# load in hg38 genomic coordinates
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene") 

# map gene symbols
mcols(exons)$gene_name <- mapIds(org.Hs.eg.db, 
                                    keys = names(exons), 
                                    column = "SYMBOL", 
                                    keytype = "ENTREZID")

# load circRNA expression data
path <- "../data/processed_cellline/common_samples/"

ciri_gcsi_sub <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F) 
ciri_gdsc_sub <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)
ciri_ccle_sub <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)

circ_gcsi_sub <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)
circ_gdsc_sub <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)
circ_ccle_sub <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F)

cfnd_gcsi_sub <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)
cfnd_gdsc_sub <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)
cfnd_ccle_sub <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F)

fcrc_gcsi_sub <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)
fcrc_gdsc_sub <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)
fcrc_ccle_sub <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F)

############################################################
# Create GRanges of circRNA genomic coordinates
############################################################

# function to create BED files
IDstoBED <- function(df) {

    circIDs <- colnames(df[,-c(1)])
    coords <- str_split(circIDs, "\\.", simplify = T) |> as.data.frame() 
    colnames(coords) <- c("chr", "start", "end")
    bed <- makeGRangesFromDataFrame(coords)
    mcols(bed)$circID <- circIDs
    
    return(bed)
}

ciri_gcsi_gr <- IDstoBED(ciri_gcsi_sub)
ciri_gdsc_gr <- IDstoBED(ciri_gdsc_sub)
ciri_ccle_gr <- IDstoBED(ciri_ccle_sub)

circ_gcsi_gr <- IDstoBED(circ_gcsi_sub)
circ_gdsc_gr <- IDstoBED(circ_gdsc_sub)
circ_ccle_gr <- IDstoBED(circ_ccle_sub)

cfnd_gcsi_gr <- IDstoBED(cfnd_gcsi_sub)
cfnd_gdsc_gr <- IDstoBED(cfnd_gdsc_sub)
cfnd_ccle_gr <- IDstoBED(cfnd_ccle_sub)

fcrc_gcsi_gr <- IDstoBED(fcrc_gcsi_sub)
fcrc_gdsc_gr <- IDstoBED(fcrc_gdsc_sub)
fcrc_ccle_gr <- IDstoBED(fcrc_ccle_sub)


############################################################
# Map circRNAs to transcript coordinates
############################################################

# function to map circRNA to genes
circToGene <- function(gr) {

    print("---Starting circToGene mapping")

    circ_genes <- list()

    for (i in 1:length(gr)) {
        genes <- c()
        circRNA = gr[i]
        overlaps <- findOverlaps(exons, circRNA)
        hits <- exons[queryHits(overlaps)]

        # ensure at least one exon is contained in circRNA transcript
        for (gene in names(hits)) {
            indiv_exons <- hits[[gene]]
            olaps <- findOverlaps(indiv_exons, circRNA, type = "within")
            if (length(olaps) > 0)  {
                genes <- c(genes, gene)
            }
        }
        circ_genes[[circRNA$circID]] <- unname(mcols(exons[genes])$gene_name)
    }

    # creating mapping df of circ to gene from circ_genes list
    mapping <- do.call(rbind, lapply(names(circ_genes), function(col) {
        if (length(circ_genes[[col]]) > 0) {
            data.frame(old_name = col, gene = circ_genes[[col]])
        } else {
            NULL
        }
    })) |> as.data.frame(stringsAsFactors = FALSE)

    return(mapping)
}


ciri_gcsi_map <- circToGene(ciri_gcsi_gr)
ciri_gdsc_map <- circToGene(ciri_gdsc_gr)
ciri_ccle_map <- circToGene(ciri_ccle_gr)

circ_gcsi_map <- circToGene(circ_gcsi_gr)
circ_gdsc_map <- circToGene(circ_gdsc_gr)
circ_ccle_map <- circToGene(circ_ccle_gr)

cfnd_gcsi_map <- circToGene(cfnd_gcsi_gr)
cfnd_gdsc_map <- circToGene(cfnd_gdsc_gr)
cfnd_ccle_map <- circToGene(cfnd_ccle_gr)

fcrc_gcsi_map <- circToGene(fcrc_gcsi_gr)
fcrc_gdsc_map <- circToGene(fcrc_gdsc_gr)
fcrc_ccle_map <- circToGene(fcrc_ccle_gr)


# checkpoint
save(ciri_gcsi_map, ciri_gdsc_map, ciri_ccle_map, 
     circ_gcsi_map, circ_gdsc_map, circ_ccle_map,
     cfnd_gcsi_map, cfnd_gdsc_map, cfnd_ccle_map,
     fcrc_gcsi_map, fcrc_gdsc_map, fcrc_ccle_map, 
     file = "../results/data/temp/circ_to_gene_map.RData")



############################################################
# Convert circRNA expression matrix labels to genes
############################################################

# function to map circRNA expression matrices to gene-level annotations
circExpToGene <- function(circ_map, circ_exp) {

    gene_exp <- circ_exp %>%
        pivot_longer(cols = -sample, names_to = "old_name", values_to = "expression") %>%
        inner_join(circ_map, by = "old_name", relationship = "many-to-many") %>%
        group_by(sample, gene) %>%
        summarize(expression = mean(expression), .groups = "drop") %>%
        pivot_wider(names_from = gene, values_from = expression) %>% 
        as.data.frame()

    return(gene_exp)

}

ciri_gcsi_ge <- circExpToGene(ciri_gcsi_map, ciri_gcsi_sub)
ciri_gdsc_ge <- circExpToGene(ciri_gdsc_map, ciri_gdsc_sub)
ciri_ccle_ge <- circExpToGene(ciri_ccle_map, ciri_ccle_sub)

circ_gcsi_ge <- circExpToGene(circ_gcsi_map, circ_gcsi_sub)
circ_gdsc_ge <- circExpToGene(circ_gdsc_map, circ_gdsc_sub)
circ_ccle_ge <- circExpToGene(circ_ccle_map, circ_ccle_sub)

cfnd_gcsi_ge <- circExpToGene(cfnd_gcsi_map, cfnd_gcsi_sub)
cfnd_gdsc_ge <- circExpToGene(cfnd_gdsc_map, cfnd_gdsc_sub)
cfnd_ccle_ge <- circExpToGene(cfnd_ccle_map, cfnd_ccle_sub)

fcrc_gcsi_ge <- circExpToGene(fcrc_gcsi_map, fcrc_gcsi_sub)
fcrc_gdsc_ge <- circExpToGene(fcrc_gdsc_map, fcrc_gdsc_sub)
fcrc_ccle_ge <- circExpToGene(fcrc_ccle_map, fcrc_ccle_sub)


# save dataframes
write.table(ciri_gcsi_ge, file = "../data/processed_cellline/GE_common_samples/CIRI2/ciri_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_gdsc_ge, file = "../data/processed_cellline/GE_common_samples/CIRI2/ciri_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_ccle_ge, file = "../data/processed_cellline/GE_common_samples/CIRI2/ciri_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(circ_gcsi_ge, file = "../data/processed_cellline/GE_common_samples/CIRCexplorer2/circ_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_gdsc_ge, file = "../data/processed_cellline/GE_common_samples/CIRCexplorer2/circ_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_ccle_ge, file = "../data/processed_cellline/GE_common_samples/CIRCexplorer2/circ_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(cfnd_gcsi_ge, file = "../data/processed_cellline/GE_common_samples/circRNA_finder/cfnd_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_gdsc_ge, file = "../data/processed_cellline/GE_common_samples/circRNA_finder/cfnd_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_ccle_ge, file = "../data/processed_cellline/GE_common_samples/circRNA_finder/cfnd_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

write.table(fcrc_gcsi_ge, file = "../data/processed_cellline/GE_common_samples/find_circ/fcrc_gcsi_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_gdsc_ge, file = "../data/processed_cellline/GE_common_samples/find_circ/fcrc_gdsc_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_ccle_ge, file = "../data/processed_cellline/GE_common_samples/find_circ/fcrc_ccle_counts.tsv", quote = F, sep = "\t", col.names = T, row.names = F)


## OLD CODE
# ========== Map circRNAs to genes ========== #

# load in annotation file from circAtlas
# link: https://ngdc.cncb.ac.cn/circatlas/links1.php - human_bed_v3.0.zip
anno <- fread("../data/human_bed_v3.0.txt")
anno$gene <- gsub("_.*", "", gsub("hsa-", "", anno$circAltas_ID))
anno$circID <- paste(anno$Chro, anno$Start, anno$End, sep = ".")


# load in annotation file from circBase
anno <- fread("../data/hsa_hg19_circRNA.bed")
anno <- anno[,c(1:3, 6)]
colnames(anno) <- c("chrom", "chromStart", "chromEnd", "strand")
anno$chromStart <- as.numeric(anno$chromStart)
anno$chromEnd <- as.numeric(anno$chromEnd)
anno <- GenomicRanges::makeGRangesFromDataFrame(anno, na.rm = TRUE)

# annotation file from circBase after LiftOver: https://genome.ucsc.edu/cgi-bin/hgLiftOver
anno <- fread("../data/hglft_genome_10a3cf_1b6e70.bed")
