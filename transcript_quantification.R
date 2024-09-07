# Create analysis-ready count matrices for circRNA, isoforms, and gene expression #

# load libraries
suppressMessages(library(data.table))
suppressMessages(library(PharmacoGx))
options(stringsAsFactors = FALSE)


###############################################
### Create cell line circRNA count matrices ###
###############################################

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


# function to match cell.id to unique.cellid from PharmacoGx
matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}


# ========== Determine biological replicates across datasets ========== #

# read in common cells
cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

#read in gcsi cell annotations
gcsi <- read.csv(file = "../data/rnaseq_meta/gcsi_rnaseq_meta.csv")
gcsi <- gcsi[,colnames(gcsi) %in% c("cellid", "alias", "Cell_line")]
gcsi$cellid <- matchToIDTable(ids=gcsi$Cell_line , tbl = cell_all, column = "GNE.cellid", returnColumn = "unique.cellid")
rownames(gcsi) <- gcsi$alias

#read in ccle cell annotations
ccle <- read.csv(file = "../data/rnaseq_meta/ccle_rnaseq_meta.csv")
ccle <- ccle[,colnames(ccle) %in% c("cellid", "Run", "Cell_Line")]
ccle$cellid <- matchToIDTable(ids=ccle$Cell_Line , tbl = cell_all, column = "CCLE.cellid", returnColumn = "unique.cellid")
rownames(ccle) <- ccle$Run

#read in gdsc cell annotations
gdsc <- read.csv("../data/rnaseq_meta/gdsc/sample_file.csv")
gdsc_mapping <- read.csv("../data/rnaseq_meta/gdsc/samples.csv")
gdsc$file_name <- gsub("\\..*", "", gdsc$file_name)
gdsc <- gdsc[gdsc$file_name %in% ciri_gdsc$V1,]
gdsc$sample_alias <- gdsc_mapping$subject_id[match(gdsc$sample_alias, gdsc_mapping$alias)]
gdsc$cellid <- matchToIDTable(ids = gdsc$sample_alias, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
rownames(gdsc) <- gdsc$file_name

# rename rownames of gcsi dataframes
ciri_gcsi$V1 <- gcsi$cellid[match(gsub("gcsi", "", ciri_gcsi$V1), rownames(gcsi))]
circ_gcsi$V1 <- gcsi$cellid[match(gsub("gcsi", "", circ_gcsi$V1), rownames(gcsi))]
cfnd_gcsi$V1 <- gcsi$cellid[match(gsub("gcsi", "", cfnd_gcsi$V1), rownames(gcsi))]
fcrc_gcsi$V1 <- gcsi$cellid[match(gsub("gcsi", "", fcrc_gcsi$V1), rownames(gcsi))]

# rename rownames of ccle dataframes
ciri_ccle$V1 <- ccle$cellid[match(ciri_ccle$V1, rownames(ccle))]
circ_ccle$V1 <- ccle$cellid[match(circ_ccle$V1, rownames(ccle))]
cfnd_ccle$V1 <- ccle$cellid[match(cfnd_ccle$V1, rownames(ccle))]
fcrc_ccle$V1 <- ccle$cellid[match(fcrc_ccle$V1, rownames(ccle))]

# rename rownames of gdsc dataframes
ciri_gdsc$V1 <- gdsc$cellid[match(ciri_gdsc$V1, rownames(gdsc))]
circ_gdsc$V1 <- gdsc$cellid[match(circ_gdsc$V1, rownames(gdsc))]
cfnd_gdsc$V1 <- gdsc$cellid[match(cfnd_gdsc$V1, rownames(gdsc))]
fcrc_gdsc$V1 <- gdsc$cellid[match(fcrc_gdsc$V1, rownames(gdsc))]


# ========== Normalize by read depth ========== #

# load in library counts
gcsi_reads <- read.table("/cluster/home/julian/circRNA_biomarker/data/raw_cellline/readCounts/gcsi.tsv")
ccle_reads <- read.table("/cluster/home/julian/circRNA_biomarker/data/raw_cellline/readCounts/ccle.tsv")
gdsc_reads <- read.table("/cluster/home/julian/circRNA_biomarker/data/raw_cellline/readCounts/gdsc.tsv")

# match sample names
gcsi_reads$sample <- gcsi$cellid[match(gcsi_reads$sample, rownames(gcsi))]
ccle_reads$sample <- ccle$cellid[match(ccle_reads$sample, rownames(ccle))]
gdsc_reads$sample <- gdsc$cellid[match(gdsc_reads$sample, rownames(gdsc))]


# function to normalize counts
cpm <- function(df, lib) {

  # save sample names
  samples <- data.frame(sample = df$V1)

  # keep only chromosomal circRNAs
  rm <- colnames(df)[-grep("chr", colnames(df))]
  if(length(rm) > 0) {df <- df[,-which(colnames(df) %in% rm)]}

  # order the samples and match order with library size
  lib <- lib[match(samples$sample, lib$sample),]

  # divide counts in each row by the corresponding library size for that sample
  df_cpm <- lapply(df, function(col) col / lib$avg_counts * 1e6)

  # re-add sample names
  df <- cbind(samples, df_cpm)

  return(df)
}

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

# function to average technical replicates               
avg_reps <- function(df) {

  # identify replicates
  reps <- df$sample[duplicated(df$sample)]
  print(reps)

  # if replicates present, loop through each replicate to average counts
  if (length(reps) > 0) {
        
    for (rep in reps) {
      print(rep)
      
      # average counts across replicates
      num_reps <- nrow(df[which(df$sample == rep),])
      tmp <- data.frame(lapply(df[which(df$sample == rep), -which(colnames(df) == "sample")], as.numeric))
      avg_counts <- c(rep, as.numeric(colSums(tmp)/num_reps)) #'x' must be numeric
      
      # replace replicates with average
      df <- df[-which(df$sample == rep),]
      df <- rbind(df, avg_counts)
    }
  }

  return(df)
        
}

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

