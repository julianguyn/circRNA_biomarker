# ========== Determine biological replicates across datasets ========== #

# read in common cells
cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

#' Function to read in pset cell annotation
#' 
#' @param filepath str. Path to cell annotation file
#' @param pset str. One of: 'gcsi', 'ccle', or 'gdsc'
#' 
load_anno <- function(filepath, pset) {
  
  anno <- read.csv(filepath)
  

}

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
gdsc[gdsc$sample_alias == "Geo",]$sample_alias = "GEO"
gdsc[gdsc$sample_alias == "NB(TU)1-10",]$sample_alias = "NB-TU-1-10"
gdsc[gdsc$sample_alias == "NTERA-2cl.D1",]$sample_alias = "NTERA-2cl-D1"
gdsc[gdsc$sample_alias == "UWB1.289",]$sample_alias = "UWB1-289"
rownames(gdsc) <- gdsc$file_name

# match sample to tissue metadata
gcsi$tissue <- gcsi_tmeta$Tissue_supergroup[match(gcsi$Cell_line, gcsi_tmeta$Cell_line)]
ccle$tissue <- str_to_title(gsub("_", " ", ccle_tmeta$tissue[match(ccle$Cell_Line, ccle_tmeta$Cell_Line)]))
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

#' CPM normalization
#' 
#' @param df data.frame. Output of circRNA_pipeline_processing.R
#' @param lib data.frame. Library sizes
#' 
get_cpm <- function(df, lib) {

    df$sample <- sub("_INPUT", "", df$sample)

    # order the samples and match order with library size
    lib <- lib[match(df$sample, lib$sample), ]
    samples <- df$sample

    # keep only chromosomal circRNAs
    rm <- colnames(df)[-grep("chr", colnames(df))]
    if(length(rm) > 0) df <- df[,-which(colnames(df) %in% rm)]

    # divide counts in each row by the corresponding library size for that sample
    df_cpm <- sapply(df, function(col) col / lib$avg_counts * 1e6) 

    # re-add sample names
    df <- data.frame(sample = samples, df_cpm)

    return(df)
}

#' Average across all technical replicates
#' 
avg_reps <- function(df) {

  # identify replicates
  reps <- df$sample[duplicated(df$sample)]

  if (length(reps) == 0) {

    print("No technical replicates")

  } else {

    print(paste("Identified", length(reps), "replicates:"))

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