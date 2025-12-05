

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

# load lung circRNA expression data
load("../data/processed_lung/circ_lung_expression.RData")

############################################################
# Format lung dataframes
############################################################

# helper function
format_lung <- function(df) {
    df <- data.frame(sample = rownames(df), df)
    rownames(df) <- NULL
    return(df)
}

ciri_polyA_gr <- format_lung(ciri_polyA)
ciri_ribo0_gr <- format_lung(ciri_ribo0)

circ_polyA_gr <- format_lung(circ_polyA)
circ_ribo0_gr <- format_lung(circ_ribo0)

cfnd_polyA_gr <- format_lung(cfnd_polyA)
cfnd_ribo0_gr <- format_lung(cfnd_ribo0)

fcrc_polyA_gr <- format_lung(fcrc_polyA)
fcrc_ribo0_gr <- format_lung(fcrc_ribo0)

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

# cell line
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

# lung
ciri_polyA_gr <- IDstoBED(ciri_polyA)
ciri_ribo0_gr <- IDstoBED(ciri_ribo0)

circ_polyA_gr <- IDstoBED(circ_polyA)
circ_ribo0_gr <- IDstoBED(circ_ribo0)

cfnd_polyA_gr <- IDstoBED(cfnd_polyA)
cfnd_ribo0_gr <- IDstoBED(cfnd_ribo0)

fcrc_polyA_gr <- IDstoBED(fcrc_polyA)
fcrc_ribo0_gr <- IDstoBED(fcrc_ribo0)


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

        # get exons in circRNA
        overlaps <- findOverlaps(unlist(exons), circRNA, type = "within")
        exons_within <- unlist(exons)[queryHits(overlaps)]
        n_exons <- length(exons_within)

        # skip if there are no exons
        if (n_exons == 0) {
            circ_genes[[circRNA$circID]] <- character(0)
            next
        }

        # table of number of exons per gene
        gene_counts <- table(names(exons_within))

        # keep genes where number of overlapping exons >= half number of exons in circRNA
        genes <- names(gene_counts[gene_counts >= floor(n_exons / 2)])
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

# cell line
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


# lung
ciri_polyA_map <- circToGene(ciri_polyA_gr)
ciri_ribo0_map <- circToGene(ciri_ribo0_gr)

circ_polyA_map <- circToGene(circ_polyA_gr)
circ_ribo0_map <- circToGene(circ_ribo0_gr)

cfnd_polyA_map <- circToGene(cfnd_polyA_gr)
cfnd_ribo0_map <- circToGene(cfnd_ribo0_gr)

fcrc_polyA_map <- circToGene(fcrc_polyA_gr)
fcrc_ribo0_map <- circToGene(fcrc_ribo0_gr)

# checkpoint
save(ciri_polyA_map, ciri_ribo0_map,
     circ_polyA_map, circ_ribo0_map,
     cfnd_polyA_map, cfnd_ribo0_map,
     fcrc_polyA_map, fcrc_ribo0_map,
     file = "../results/data/temp/circ_lung_to_gene_map.RData")

############################################################
# Convert circRNA expression matrix labels to genes
############################################################

# function to map circRNA expression matrices to gene-level annotations
circExpToGene <- function(circ_map, circ_exp, dir, label) {

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
