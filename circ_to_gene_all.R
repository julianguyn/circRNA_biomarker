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
path <- "../data/processed_cellline/all_samples/"

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
# Identify cell lines to keep
############################################################

# find intersecting cell lines
gcsi_ccle <- intersect(ciri_gcsi_sub$sample,ciri_ccle_sub$sample)   # 481 samples
gcsi_gdsc <- intersect(ciri_gcsi_sub$sample,ciri_gdsc_sub$sample)   # 112 samples
ccle_gdsc <- intersect(ciri_ccle_sub$sample,ciri_gdsc_sub$sample)   # 91 samples

to_keep <- unique(c(gcsi_ccle, gcsi_gdsc, ccle_gdsc))               # 588 total
save(to_keep, file = "../data/temp/all_intersect_cells.RData")

############################################################
# Filter samples
############################################################

# function to filter circ expressiond dataframes 
filter_circ <- function(circ_counts) {
    circ_counts <- circ_counts[circ_counts$sample %in% to_keep,]
    return(circ_counts)
}

ciri_gcsi_sub <- filter_circ(ciri_gcsi_sub)
ciri_gdsc_sub <- filter_circ(ciri_gdsc_sub)
ciri_ccle_sub <- filter_circ(ciri_ccle_sub)

circ_gcsi_sub <- filter_circ(circ_gcsi_sub)
circ_gdsc_sub <- filter_circ(circ_gdsc_sub)
circ_ccle_sub <- filter_circ(circ_ccle_sub)

cfnd_gcsi_sub <- filter_circ(cfnd_gcsi_sub)
cfnd_gdsc_sub <- filter_circ(cfnd_gdsc_sub)
cfnd_ccle_sub <- filter_circ(cfnd_ccle_sub)

fcrc_gcsi_sub <- filter_circ(fcrc_gcsi_sub)
fcrc_gdsc_sub <- filter_circ(fcrc_gdsc_sub)
fcrc_ccle_sub <- filter_circ(fcrc_ccle_sub)


############################################################
# Create GRanges of circRNA genomic coordinates
############################################################

# function to create BED files
IDstoBED <- function(df) {

    circIDs <- colnames(df[,-c(1)])
    coords <- str_split(circIDs, "\\.", simplify = T) |> as.data.frame() 
    colnames(coords) <- c("chr", "start", "end")
    coords <- coords[as.numeric(coords$end) >= as.numeric(coords$start) - 1, ]
    bed <- makeGRangesFromDataFrame(coords)
    mcols(bed)$circID <- paste(coords$chr, coords$start, coords$end, sep = ".")
    
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
path = "../data/processed_cellline/GE_all_samples/"

write.table(ciri_gcsi_ge, file = paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_gdsc_ge, file = paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(ciri_ccle_ge, file = paste0(path, "CIRI2/ciri_ccle_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

write.table(circ_gcsi_ge, file = paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_gdsc_ge, file = paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(circ_ccle_ge, file = paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

write.table(cfnd_gcsi_ge, file = paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_gdsc_ge, file = paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(cfnd_ccle_ge, file = paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

write.table(fcrc_gcsi_ge, file = paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_gdsc_ge, file = paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
write.table(fcrc_ccle_ge, file = paste0(path, "find_circ/fcrc_ccle_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

