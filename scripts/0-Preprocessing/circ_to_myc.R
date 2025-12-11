# load libraries
suppressPackageStartupMessages({
    library(data.table)
})

############################################################
# Load in circRNA expression data
############################################################

path <- "../data/processed_cellline/all_samples/"

# helper function to load in expression matrices
load_counts <- function(path, filename) {
    counts <- fread(paste0(path, filename), data.table = FALSE)
    return(counts)
}

# load circRNA expression data
ciri_gcsi <- load_counts(path, "CIRI2/ciri_gcsi_counts.tsv")
ciri_gdsc <- load_counts(path, "CIRI2/ciri_gdsc_counts.tsv")
ciri_ccle <- load_counts(path, "CIRI2/ciri_ccle_counts.tsv")

circ_gcsi <- load_counts(path, "CIRCexplorer2/circ_gcsi_counts.tsv")
circ_gdsc <- load_counts(path, "CIRCexplorer2/circ_gdsc_counts.tsv")
circ_ccle <- load_counts(path, "CIRCexplorer2/circ_ccle_counts.tsv")

cfnd_gcsi <- load_counts(path, "circRNA_finder/cfnd_gcsi_counts.tsv")
cfnd_gdsc <- load_counts(path, "circRNA_finder/cfnd_gdsc_counts.tsv")
cfnd_ccle <- load_counts(path, "circRNA_finder/cfnd_ccle_counts.tsv")

fcrc_gcsi <- load_counts(path, "find_circ/fcrc_gcsi_counts.tsv")
fcrc_gdsc <- load_counts(path, "find_circ/fcrc_gdsc_counts.tsv")
fcrc_ccle <- load_counts(path, "find_circ/fcrc_ccle_counts.tsv")


############################################################
# Get all circRNAs on chr8
############################################################

all_circ <- c(
    colnames(ciri_gcsi),
    colnames(ciri_gdsc),
    colnames(ciri_ccle),
    colnames(circ_gcsi),
    colnames(circ_gdsc),
    colnames(circ_ccle),
    colnames(cfnd_gcsi),
    colnames(cfnd_gdsc),
    colnames(cfnd_ccle),
    colnames(fcrc_gcsi),
    colnames(fcrc_gdsc),
    colnames(fcrc_ccle)
)

# get all circRNAs on chr8
chr8 <- gsub("chr8\\.", "", all_circ[grep("chr8", all_circ)])
split_vals <- strsplit(chr8, "\\.")
df <- data.frame(
    circ = paste0("chr8.", chr8),
    start = as.numeric(sapply(split_vals, `[`, 1)),
    end = as.numeric(sapply(split_vals, `[`, 2)),
    myc = 0
)

############################################################
# Identify circRNA mapping to MYC
############################################################

# get MYC bases (MYC gene: 8:127735434-127742951)
myc_start <- 127735434
myc_end <- 127742951
myc <- myc_start:myc_end

# save number of overlapping bases per circRNA on chr8
for (i in seq_along(df$circ)) {
    circ <- df$start[i]:df$end[i]
    if (sum(circ %in% myc) > 0) {
        df$myc[i] <- sum(circ %in% myc)
    }
}

# identify circRNAs that overlapp 100% with the MYC gene
df$total <- df$end - df$start + 1
df$prop <- df$myc / df$total * 100

# save all circMYCs
df <- df[df$prop == 100,] # n = 508
save(df, file = "../results/data/circ_mapped_to_myc.RData")


############################################################
# Keep circRNA overlapping 100% to MYC
############################################################

# helper function to filter circ exp dataframes
filter_myc <- function(circ) {

    rownames(circ) <- circ$sample
    circ <- circ[,colnames(circ) %in% df$circ]
    print(ncol(circ))
    circ <- data.frame(MYC = rowMeans(circ))

    return(circ)
}

ciri_gcsi <- filter_myc(ciri_gcsi) #8
ciri_gdsc <- filter_myc(ciri_gdsc) #2
ciri_ccle <- filter_myc(ciri_ccle) #9

circ_gcsi <- filter_myc(circ_gcsi) #85
circ_gdsc <- filter_myc(circ_gdsc) #3
circ_ccle <- filter_myc(circ_ccle) #19

cfnd_gcsi <- filter_myc(cfnd_gcsi) #116
cfnd_gdsc <- filter_myc(cfnd_gdsc) #3
cfnd_ccle <- filter_myc(cfnd_ccle) #155

fcrc_gcsi <- filter_myc(fcrc_gcsi) #41
fcrc_gdsc <- filter_myc(fcrc_gdsc) #11
fcrc_ccle <- filter_myc(fcrc_ccle) #56

# save dataframes
save(ciri_gcsi, ciri_gdsc, ciri_ccle,
     circ_gcsi, circ_gdsc, circ_ccle,
     cfnd_gcsi, cfnd_gdsc, cfnd_ccle,
     fcrc_gcsi, fcrc_gdsc, fcrc_ccle,
     file = "../results/data/circ_myc.RData")