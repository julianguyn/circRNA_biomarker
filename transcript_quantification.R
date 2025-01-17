# Create analysis-ready count matrices for circRNA, isoforms, and gene expression #

# load libraries
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggvenn))
suppressMessages(library(PharmacoGx))
suppressMessages(library(stringr))
suppressMessages(library(ComplexHeatmap))

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


# load in metadata for tissue-specific analysis
#gcsi <- readRDS("gCSI.rds") 
#gcsi <- updateObject(gcsi)
#gcsi_metadata <- gcsi@molecularProfiles$Kallisto_0.46.1.rnaseq@colData #extract metadata of cell lines
#write.table(gcsi_metadata, file = "../data/rnaseq_meta/tissue_meta/gcsi_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#gdsc <- readRDS("GDSC2-8.2.rds")
#gdsc <- updateObject(gdsc)
#gdsc_metadata <- gdsc@molecularProfiles$Kallisto_0.46.1.rnaseq@colData
#write.table(gdsc_metadata, file = "../data/rnaseq_meta/tissue_meta/gdsc_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

#ccle <- readRDS("CCLE.rds")
#ccle <- updateObject(ccle)
#ccle_metadata <- ccle@molecularProfiles$Kallisto_0.46.1.rnaseq@colData
#write.table(ccle_metadata, file = "../data/rnaseq_meta/tissue_meta/ccle_metadata.tsv", quote = F, sep = "\t", col.names = T, row.names = F)

# load in tissue metadata
gcsi_tmeta <- fread("../data/rnaseq_meta/tissue_meta/gcsi_metadata.tsv")
ccle_tmeta <- fread("../data/rnaseq_meta/tissue_meta/ccle_metadata.tsv")
gdsc_tmeta <- fread("../data/rnaseq_meta/tissue_meta/gdsc_metadata.tsv")


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

# match gdsc annotations
gdsc$cellid <- matchToIDTable(ids = gdsc$sample_alias, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
gdsc[is.na(gdsc$cellid),]$cellid <- matchToIDTable(ids = gdsc[is.na(gdsc$cellid),]$sample_alias, tbl = cell_all, column = "GDSC1000.cellid", returnColumn = "unique.cellid")
gdsc[gdsc$sample_alias == "NCI-H820",]$cellid = "NCI-H820"
gdsc[gdsc$sample_alias == "PC-3-JPC-3",]$cellid = "PC-3 [Human lung carcinoma]"
gdsc[gdsc$sample_alias == "Geo",]$cellid = "GEO"
gdsc[gdsc$sample_alias == "NTERA-2cl.D1",]$cellid = "NTERA-2"
rownames(gdsc) <- gdsc$file_name

# match sample to tissue metadata
gcsi$tissue <- gcsi_tmeta$Tissue_supergroup[match(gcsi$Cell_line, gcsi_tmeta$Cell_line)]
ccle$tissue <- str_to_title(gsub("_", " ", ccle_tmeta$tissue[match(ccle$Cell_Line, ccle_tmeta$Cell_Line)]))
gdsc[gdsc$sample_alias == "Geo",]$sample_alias = "GEO"
gdsc[gdsc$sample_alias == "NB(TU)1-10",]$sample_alias = "NB-TU-1-10"
gdsc[gdsc$sample_alias == "NTERA-2cl.D1",]$sample_alias = "NTERA-2cl-D1"
gdsc[gdsc$sample_alias == "UWB1.289",]$sample_alias = "UWB1-289"
gdsc$tissue <- str_to_title(gdsc_tmeta$Factor.Value.organism.part.[match(gdsc$sample_alias, gdsc_tmeta$Source.Name)])

# save metadata file
save(gcsi, ccle, gdsc, file = "../results/data/tissue-metadata.RData")

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



####################################################
### Figure 1: circRNA Quantification Comparisons ###
####################################################

## TODO: do this for common samples

# set palette for plotting
pal = c("#51C7AD", "#392C57", "#3670A0")

# remove sample names column for downstream analysis
ciri_gcsi <- ciri_gcsi[,-which(colnames(ciri_gcsi) %in% c("sample"))]
ciri_ccle <- ciri_ccle[,-which(colnames(ciri_ccle) %in% c("sample"))]
ciri_gdsc <- ciri_gdsc[,-which(colnames(ciri_gdsc) %in% c("sample"))]

circ_gcsi <- circ_gcsi[,-which(colnames(circ_gcsi) %in% c("sample"))]
circ_ccle <- circ_ccle[,-which(colnames(circ_ccle) %in% c("sample"))]
circ_gdsc <- circ_gdsc[,-which(colnames(circ_gdsc) %in% c("sample"))]

cfnd_gcsi <- cfnd_gcsi[,-which(colnames(cfnd_gcsi) %in% c("sample"))]
cfnd_ccle <- cfnd_ccle[,-which(colnames(cfnd_ccle) %in% c("sample"))]
cfnd_gdsc <- cfnd_gdsc[,-which(colnames(cfnd_gdsc) %in% c("sample"))]

fcrc_gcsi <- fcrc_gcsi[,-which(colnames(fcrc_gcsi) %in% c("sample"))]
fcrc_ccle <- fcrc_ccle[,-which(colnames(fcrc_ccle) %in% c("sample"))]
fcrc_gdsc <- fcrc_gdsc[,-which(colnames(fcrc_gdsc) %in% c("sample"))]

# ========== Unique Transcript Detection: Dataset Comparison per Pipeline ========== #

# create list object of transcripts
ciri_transcripts <- list(gCSI = colnames(ciri_gcsi), CCLE = colnames(ciri_ccle), GDSC2 = colnames(ciri_gdsc))
circ_transcripts <- list(gCSI = colnames(circ_gcsi), CCLE = colnames(circ_ccle), GDSC2 = colnames(circ_gdsc))
cfnd_transcripts <- list(gCSI = colnames(cfnd_gcsi), CCLE = colnames(cfnd_ccle), GDSC2 = colnames(cfnd_gdsc))
fcrc_transcripts <- list(gCSI = colnames(fcrc_gcsi), CCLE = colnames(fcrc_ccle), GDSC2 = colnames(fcrc_gdsc))

# plot venn diagram
p1 <- ggvenn(ciri_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "CIRI2")
p2 <- ggvenn(circ_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "CIRCexplorer2")
p3 <- ggvenn(cfnd_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "circRNA_finder")
p4 <- ggvenn(fcrc_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "find_circ")

png("../results/figures/figure1/venndiagram_per_pipeline.png", width=500, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = FALSE)
dev.off()

png("../results/figures/figure1/common_venndiagram_per_pipeline.png", width=500, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = FALSE)
dev.off()

# ========== Unique Transcript Detection: Pipeline Comparison ========== #

# create list object of transcripts for venn diagram
all_comparison <- list(CIRI2 = c(colnames(ciri_gcsi), colnames(ciri_ccle), colnames(ciri_gdsc)), 
                       CIRCexplorer2 = c(colnames(circ_gcsi), colnames(circ_ccle), colnames(circ_gdsc)), 
                       circRNA_finder = c(colnames(cfnd_gcsi), colnames(cfnd_ccle), colnames(cfnd_gdsc)),
                       find_circ = c(colnames(fcrc_gcsi), colnames(fcrc_ccle), colnames(fcrc_gdsc)))

# plot venn diagram
png("../results/figures/figure1/venndiagram.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggvenn(all_comparison, 
        fill_color = c("#E8F6B1", "#A5DBB7", "#2088BC", "#26479D"),
        stroke_size = 0.5, set_name_size = 4)
dev.off()


# create list object of transcripts for upset plot
set.seed(123)
toPlot <- make_comb_mat(list(
            gCSI_CIRI2 = colnames(ciri_gcsi), CCLE_CIRI2 = colnames(ciri_ccle), GDSC_CIRI2 = colnames(ciri_gdsc),
            gCSI_CIRCexplorer2 = colnames(circ_gcsi), CCLE_CIRCexplorer2 = colnames(circ_ccle), GDSC_CIRCexplorer2 = colnames(circ_gdsc),
            gCSI_circRNA_finder = colnames(cfnd_gcsi), CCLE_circRNA_finder = colnames(cfnd_ccle), GDSC_circRNA_finder = colnames(cfnd_gdsc),
            gCSI_find_circ = colnames(fcrc_gcsi), CCLE_find_circ = colnames(fcrc_ccle), GDSC_find_circ = colnames(fcrc_gdsc)))

# remove combinations of less than 100 pairs
toPlot <- toPlot[comb_size(toPlot) >= 100]

# upset plot
pdf("../results/figures/figure1/upset_by_pipeline.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "CCLE_CIRI2", "GDSC_CIRI2", "gCSI_CIRCexplorer2", "CCLE_CIRCexplorer2", "GDSC_CIRCexplorer2",
                            "gCSI_circRNA_finder", "CCLE_circRNA_finder", "GDSC_circRNA_finder", "gCSI_find_circ", "CCLE_find_circ", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()

pdf("../results/figures/figure1/upset_by_pset.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "gCSI_CIRCexplorer2", "gCSI_circRNA_finder", "gCSI_find_circ",
                            "CCLE_CIRCexplorer2", "CCLE_CIRI2", "CCLE_circRNA_finder", "CCLE_find_circ",
                            "GDSC_CIRI2", "GDSC_CIRCexplorer2", "GDSC_circRNA_finder", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()

pdf("../results/figures/figure1/common_upset_by_pipeline.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "CCLE_CIRI2", "GDSC_CIRI2", "gCSI_CIRCexplorer2", "CCLE_CIRCexplorer2", "GDSC_CIRCexplorer2",
                            "gCSI_circRNA_finder", "CCLE_circRNA_finder", "GDSC_circRNA_finder", "gCSI_find_circ", "CCLE_find_circ", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()

pdf("../results/figures/figure1/common_upset_by_pset.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "gCSI_CIRCexplorer2", "gCSI_circRNA_finder", "gCSI_find_circ",
                            "CCLE_CIRCexplorer2", "CCLE_CIRI2", "CCLE_circRNA_finder", "CCLE_find_circ",
                            "GDSC_CIRI2", "GDSC_CIRCexplorer2", "GDSC_circRNA_finder", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()



# ========== Transcript Quantification: Dataset Comparison per Pipeline ========== #


# TODO: remove
ciri_gdsc[is.na(ciri_gdsc)] <- 0
circ_gdsc[is.na(circ_gdsc)] <- 0
cfnd_gdsc[is.na(cfnd_gdsc)] <- 0
fcrc_gdsc[is.na(fcrc_gdsc)] <- 0

# create data frame of counts for plotting
df <- data.frame(Count = c(sum(ciri_gcsi), sum(ciri_ccle), sum(ciri_gdsc),
                           sum(circ_gcsi), sum(circ_ccle), sum(circ_gdsc), 
                           sum(cfnd_gcsi), sum(cfnd_ccle), sum(cfnd_gdsc), 
                           sum(fcrc_gcsi), sum(fcrc_ccle), sum(fcrc_gdsc)),
                PSet = c(rep(c("gCSI", "CCLE", "GDSC2"), 4)),
                Pipeline = c(rep("CIRI2", 3), rep("CIRCexplorer2", 3), rep("circRNA_finder", 3), rep("find_circ", 3)))
df$Pipeline <- factor(df$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
df$PSet <- factor(df$PSet, levels = c("gCSI", "CCLE", "GDSC2"))

# plot bar plot of counts
png("../results/figures/figure1/counts.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(df, aes(x = Pipeline, y = Count, fill = PSet)) + geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_fill_manual(values=pal, limits=c("gCSI", "CCLE", "GDSC2")) + 
  scale_y_continuous(limits = c(0, 240000), expand=c(0,0))  + theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_counts.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(df, aes(x = Pipeline, y = log2(Count), fill = PSet)) + geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_fill_manual(values=pal, limits=c("gCSI", "CCLE", "GDSC2")) + 
  scale_y_continuous(limits = c(0, 13), expand=c(0,0))  + theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.4, 'cm')) + labs(y = "Log2 Normalized Counts")
dev.off()


# ========== Distribution of Transcript Expression ========== #

suppressWarnings(ciri_gcsi <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gcsi_counts.tsv", data.table = F))
suppressWarnings(ciri_gdsc <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gdsc_counts.tsv", data.table = F))
suppressWarnings(ciri_ccle <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_ccle_counts.tsv", data.table = F))

suppressWarnings(circ_gcsi <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F))
suppressWarnings(circ_gdsc <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F))
suppressWarnings(circ_ccle <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_ccle_counts.tsv", data.table = F))

suppressWarnings(cfnd_gcsi <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F))
suppressWarnings(cfnd_gdsc <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F))
suppressWarnings(cfnd_ccle <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_ccle_counts.tsv", data.table = F))

suppressWarnings(fcrc_gcsi <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gcsi_counts.tsv", data.table = F))
suppressWarnings(fcrc_gdsc <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gdsc_counts.tsv", data.table = F))
suppressWarnings(fcrc_ccle <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_ccle_counts.tsv", data.table = F))


# remove sample names column for downstream analysis
ciri_gcsi <- ciri_gcsi[,-which(colnames(ciri_gcsi) %in% c("sample"))]
ciri_ccle <- ciri_ccle[,-which(colnames(ciri_ccle) %in% c("sample"))]
ciri_gdsc <- ciri_gdsc[,-which(colnames(ciri_gdsc) %in% c("sample"))]

circ_gcsi <- circ_gcsi[,-which(colnames(circ_gcsi) %in% c("sample"))]
circ_ccle <- circ_ccle[,-which(colnames(circ_ccle) %in% c("sample"))]
circ_gdsc <- circ_gdsc[,-which(colnames(circ_gdsc) %in% c("sample"))]

cfnd_gcsi <- cfnd_gcsi[,-which(colnames(cfnd_gcsi) %in% c("sample"))]
cfnd_ccle <- cfnd_ccle[,-which(colnames(cfnd_ccle) %in% c("sample"))]
cfnd_gdsc <- cfnd_gdsc[,-which(colnames(cfnd_gdsc) %in% c("sample"))]

fcrc_gcsi <- fcrc_gcsi[,-which(colnames(fcrc_gcsi) %in% c("sample"))]
fcrc_ccle <- fcrc_ccle[,-which(colnames(fcrc_ccle) %in% c("sample"))]
fcrc_gdsc <- fcrc_gdsc[,-which(colnames(fcrc_gdsc) %in% c("sample"))]


# create data frame of expression values
print("making expression")
Expression = c(as.vector(unlist(ciri_gcsi)), as.vector(unlist(ciri_ccle)), as.vector(unlist(ciri_gdsc)),
                                    as.vector(unlist(circ_gcsi)), as.vector(unlist(circ_ccle)), as.vector(unlist(circ_gdsc)),
                                    as.vector(unlist(cfnd_gcsi)), as.vector(unlist(cfnd_ccle)), as.vector(unlist(cfnd_gdsc)),
                                    as.vector(unlist(fcrc_gcsi)), as.vector(unlist(fcrc_ccle)), as.vector(unlist(fcrc_gdsc)))
print("making label")
Label = c(rep("gCSI-CIRI2", prod(dim(ciri_gcsi))), rep("CCLE-CIRI2", prod(dim(ciri_ccle))), rep("GDSC2-CIRI2", prod(dim(ciri_gdsc))),
                               rep("gCSI-CIRCexplorer2", prod(dim(circ_gcsi))), rep("CCLE-CIRCexplorer2", prod(dim(circ_ccle))), rep("GDSC2-CIRCexplorer2", prod(dim(circ_gdsc))),
                               rep("gCSI-circRNA_finder", prod(dim(cfnd_gcsi))), rep("CCLE-circRNA_finder", prod(dim(cfnd_ccle))), rep("GDSC2-circRNA_finder", prod(dim(cfnd_gdsc))),
                               rep("gCSI-find_circ", prod(dim(fcrc_gcsi))), rep("CCLE-find_circ", prod(dim(fcrc_ccle))), rep("GDSC2-find_circ", prod(dim(fcrc_gdsc))))
                
toPlot <- data.frame(Expression, Label)
dim(toPlot)
toPlot$Pipeline <- gsub(".*-", "", toPlot$Label)
toPlot$PSet <- gsub("-.*", "", toPlot$Label)
toPlot$Pipeline <- factor(toPlot$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
toPlot$PSet <- factor(toPlot$PSet, levels = c("gCSI", "CCLE", "GDSC2"))


max = max(toPlot$Expression)

png("../results/figures/figure1/count_dist_all.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pipeline.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

# COMMON:
png("../results/figures/figure1/common_count_dist_all.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pipeline.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

# compute distribution of 0 values in normalized circRNA counts
library(dplyr)
toPlot$Group <- toPlot$Expression == 0
toPlot_counts <- toPlot %>% group_by(Pipeline, Group) %>% summarise(Count = n()) 
toPlot_counts

toPlot <- toPlot[toPlot$Expression != 0,]
save(toPlot, toPlot_counts, file = "temp.RData")


toPlot_counts <- as.data.frame(toPlot_counts)
# pie chart
p1 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "CIRI2",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("CIRI2")
p2 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "CIRCexplorer2",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("CIRCexplorer2")
p3 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "circRNA_finder",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("circRNA_finder")
p4 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "find_circ",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("find_circ")
png("count_zero_pipeline.png", width=100, height=100, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()

png("../results/figures/figure1/common_count_zero_pipeline.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot_counts, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    facet_grid(.~Pipeline) + theme(legend.key.size = unit(0.4, 'cm')) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")
dev.off()


toPlot_counts <- toPlot %>% group_by(PSet, Group) %>% summarise(Count = n()) 
toPlot_counts

# pie chart
png("../results/figures/figure1/count_zero_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot_counts, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    facet_grid(.~PSet) + theme(legend.key.size = unit(0.4, 'cm')) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")
dev.off()

png("../results/figures/figure1/common_count_zero_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot_counts, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    facet_grid(.~PSet) + theme(legend.key.size = unit(0.4, 'cm')) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")
dev.off()


# remove all 0s
toPlot <- toPlot[toPlot$Expression != 0,]

max = max(toPlot$Expression)

png("count_dist_common_nozero.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Pipeline, y = log2(Expression))) + 
    geom_violin(aes(fill = Pipeline), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.3) +
    theme_classic() + labs(x = "", fill = "", y = "log2 Normalized circRNA Count Value") +
    scale_fill_manual(values = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    coord_flip()
dev.off()

ggplot(toPlot, aes(x = log2(Expression))) + geom_density(aes(fill = Pipeline), alpha = 0.5, size = 1) + 
        theme_classic() + facet_grid(factor(Pipeline)~.) +
        scale_fill_manual(values = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "log2 Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))

png("count_dist_common_nozero.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = log2(Expression))) + geom_density(aes(fill = Pipeline), alpha = 0.5) + 
        theme_classic() + 
        scale_fill_manual(values = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "log2 Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_all_nozero.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = log2(Expression))) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "log2 Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pipeline_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pset_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()


#COMMON:
png("../results/figures/figure1/common_count_dist_all_nozero.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pipeline_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pset_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()


##TODO: Add code for lung, isoform, and gene expression processing, normalization