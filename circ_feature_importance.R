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
    library(caret)
})

set.seed(200) 

############################################################
# Load in data 
############################################################

# load circRNA expression data (subsetted in si_distribution.R)
load("../results/data/temp/circ_stability_subsetdf.RData")


############################################################
# Initiate feature dataframe
############################################################

# function to compute median expression and length
getFeats <- function(df) {

    # initiate dataframe to store results
    res <- data.frame(matrix(nrow = ncol(df), ncol = 4))
    rownames(res) <- colnames(df)
    colnames(res) <- c("MedianExp", "Length", "NExons", "GC")

    # compute median expression
    res$MedianExp <- colSums(df)/48 |> as.numeric()

    # compute length
    coords <- str_split(rownames(res), "\\.", simplify = T) |> as.data.frame()
    res$Length <- as.numeric(coords$V3) - as.numeric(coords$V2) + 1

    return(res)
}

ciri_gcsi_ft <- getFeats(ciri_gcsi)
ciri_ccle_ft <- getFeats(ciri_ccle)
ciri_gdsc_ft <- getFeats(ciri_gdsc)

circ_gcsi_ft <- getFeats(circ_gcsi)
circ_ccle_ft <- getFeats(circ_ccle)
circ_gdsc_ft <- getFeats(circ_gdsc)

cfnd_gcsi_ft <- getFeats(cfnd_gcsi)
cfnd_ccle_ft <- getFeats(cfnd_ccle)
cfnd_gdsc_ft <- getFeats(cfnd_gdsc)

fcrc_gcsi_ft <- getFeats(fcrc_gcsi)
fcrc_ccle_ft <- getFeats(fcrc_ccle)
fcrc_gdsc_ft <- getFeats(fcrc_gdsc)


############################################################
# Extract GC%
############################################################

# load in GRCh38 reference genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# function to compute gc content
computeGC <- function(res) {

    coords <- str_split(rownames(res), "\\.", simplify = T) |> as.data.frame() 
    colnames(coords) <- c("chr", "start", "end")
    for (i in 1:nrow(coords)) {
        transcript_seq <- BSgenome::getSeq(genome, names = coords$chr[i], start = as.numeric(coords$start[i]), end = as.numeric(coords$end[i]))
        gc_percent <- letterFrequency(transcript_seq, letters = c("G", "C"), as.prob = TRUE) |> sum() * 100
        res$GC[i] <- gc_percent
    }
    return(res)
}

ciri_gcsi_ft <- computeGC(ciri_gcsi_ft)
ciri_ccle_ft <- computeGC(ciri_ccle_ft)
ciri_gdsc_ft <- computeGC(ciri_gdsc_ft)

circ_gcsi_ft <- computeGC(circ_gcsi_ft)
circ_ccle_ft <- computeGC(circ_ccle_ft)
circ_gdsc_ft <- computeGC(circ_gdsc_ft)

cfnd_gcsi_ft <- computeGC(cfnd_gcsi_ft)
cfnd_ccle_ft <- computeGC(cfnd_ccle_ft)
cfnd_gdsc_ft <- computeGC(cfnd_gdsc_ft)

fcrc_gcsi_ft <- computeGC(fcrc_gcsi_ft)
fcrc_ccle_ft <- computeGC(fcrc_ccle_ft)
fcrc_gdsc_ft <- computeGC(fcrc_gdsc_ft)


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

ciri_gcsi_gr <- IDstoBED(ciri_gcsi)
ciri_gdsc_gr <- IDstoBED(ciri_gdsc)
ciri_ccle_gr <- IDstoBED(ciri_ccle)

circ_gcsi_gr <- IDstoBED(circ_gcsi)
circ_gdsc_gr <- IDstoBED(circ_gdsc)
circ_ccle_gr <- IDstoBED(circ_ccle)

cfnd_gcsi_gr <- IDstoBED(cfnd_gcsi)
cfnd_gdsc_gr <- IDstoBED(cfnd_gdsc)
cfnd_ccle_gr <- IDstoBED(cfnd_ccle)

fcrc_gcsi_gr <- IDstoBED(fcrc_gcsi)
fcrc_gdsc_gr <- IDstoBED(fcrc_gdsc)
fcrc_ccle_gr <- IDstoBED(fcrc_ccle)


############################################################
# Extract number of exons
############################################################

# load in hg38 genomic coordinates
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene") 

# map gene symbols
mcols(exons)$gene_name <- mapIds(org.Hs.eg.db, 
                                    keys = names(exons), 
                                    column = "SYMBOL", 
                                    keytype = "ENTREZID")

# function to map circRNA to genes
circToGene <- function(gr, res) {

    print("---Starting circToGene mapping")

    circ_genes <- list()

    for (i in seq_along(gr)) {
        circRNA = gr[i]
        overlaps <- findOverlaps(exons, circRNA)
        hits <- exons[queryHits(overlaps)]
        unique_exons <- GRanges()

        # ensure at least one exon is contained in circRNA transcript
        for (gene in names(hits)) {
            indiv_exons <- hits[[gene]]
            olaps <- findOverlaps(indiv_exons, circRNA, type = "within")
            if (length(olaps) > 0)  {
                unique_exons <- c(unique_exons, indiv_exons[queryHits(olaps)])
            }
        }
        unique_exons <- reduce(unique_exons)
        res$NExons[i] <- length(unique_exons)
    }

    return(res)
}

ciri_gcsi_ft <- circToGene(ciri_gcsi_gr, ciri_gcsi_ft)
ciri_gdsc_ft <- circToGene(ciri_gdsc_gr, ciri_gdsc_ft)
ciri_ccle_ft <- circToGene(ciri_ccle_gr, ciri_ccle_ft)

circ_gcsi_ft <- circToGene(circ_gcsi_gr, circ_gcsi_ft)
circ_gdsc_ft <- circToGene(circ_gdsc_gr, circ_gdsc_ft)
circ_ccle_ft <- circToGene(circ_ccle_gr, circ_ccle_ft)

cfnd_gcsi_ft <- circToGene(cfnd_gcsi_gr, cfnd_gcsi_ft)
cfnd_gdsc_ft <- circToGene(cfnd_gdsc_gr, cfnd_gdsc_ft)
cfnd_ccle_ft <- circToGene(cfnd_ccle_gr, cfnd_ccle_ft)

fcrc_gcsi_ft <- circToGene(fcrc_gcsi_gr, fcrc_gcsi_ft)
fcrc_gdsc_ft <- circToGene(fcrc_gdsc_gr, fcrc_gdsc_ft)
fcrc_ccle_ft <- circToGene(fcrc_ccle_gr, fcrc_ccle_ft)


############################################################
# Save results
############################################################

save(ciri_gcsi_ft, ciri_ccle_ft, ciri_gdsc_ft, 
     circ_gcsi_ft, circ_ccle_ft, circ_gdsc_ft,
     cfnd_gcsi_ft, cfnd_ccle_ft, cfnd_gdsc_ft, 
     fcrc_gcsi_ft, fcrc_ccle_ft, fcrc_gdsc_ft,
     file = "../results/data/circ_stability_features.RData")


# Miscellaneous
# create TxDB from the gencode annotation
#gtf_file <- "GENCODE.v45.annotation.gtf"
#txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
#saveDb(txdb, "TxDb.GENCODE.v45.sqlite")

#txdb <- loadDb("TxDb.GENCODE.v45.sqlite")
#exons <- exonsBy(txdb, by = "gene")


############################################################
# Load in circRNA stability and features
############################################################

load("../results/data/temp/circ_stability.RData")
load("../results/data/circ_stability_features.RData")