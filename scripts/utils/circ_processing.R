#' Helper function to match cell.id to unique.cellid from PharmacoGx
#' 
matchToIDTable <- function(ids, tbl, column, returnColumn = "unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}

#' Function to read in pset cell annotation
#' 
#' @param filepath str. Path to cell annotation file
#' @param pset str. One of: 'gcsi', 'ccle', or 'gdsc'
#' 
load_anno <- function(filepath, pset) {

  # read in cell annotation
  cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
  name <- switch(pset, gcsi = "alias", ccle = "Run", gdsc = "file_name")
  
  anno <- read.csv(filepath)
  if (pset %in% c("gcsi", "ccle")) {
    cell <- switch(pset, gcsi = "Cell_line", ccle = "Cell_Line")
    col <- switch(pset, gcsi = "GNE.cellid", ccle = "CCLE.cellid")
    anno <- anno[,colnames(anno) %in% c("cellid", name, "Cell_Line", "Cell_line")]
    anno$cellid <- matchToIDTable(ids=anno[[cell]], tbl = cell_all, column = col, returnColumn = "unique.cellid")
  } else {
    anno$file_name <- gsub("\\..*", "", anno$file_name)
    mapping <- read.csv("../data/rnaseq_meta/gdsc/samples.csv")
    anno$sample_alias <- mapping$subject_id[match(anno$sample_alias, mapping$alias)]
    anno$cellid <- matchToIDTable(ids = anno$sample_alias, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
    anno[is.na(anno$cellid),]$cellid <- matchToIDTable(ids = anno[is.na(anno$cellid),]$sample_alias, tbl = cell_all, column = "GDSC1000.cellid", returnColumn = "unique.cellid")
    # manual mapping for gdsc
    anno[anno$sample_alias == "NCI-H820",]$cellid <- "NCI-H820"
    anno[anno$sample_alias == "PC-3-JPC-3",]$cellid <- "PC-3 [Human lung carcinoma]"
    anno[anno$sample_alias == "Geo",]$cellid <- "GEO"
    anno[anno$sample_alias == "NTERA-2cl.D1",]$cellid <- "NTERA-2"
    anno[anno$sample_alias == "Geo",]$sample_alias <- "GEO"
    anno[anno$sample_alias == "NB(TU)1-10",]$sample_alias <- "NB-TU-1-10"
    anno[anno$sample_alias == "NTERA-2cl.D1",]$sample_alias <- "NTERA-2cl-D1"
    anno[anno$sample_alias == "UWB1.289",]$sample_alias <- "UWB1-289"
  }
  rownames(anno) <- anno[[name]]

  # match sample to tissue metadata
  tissue_meta <- fread(paste0("../data/rnaseq_meta/tissue_meta/", pset, "_metadata.tsv"))

  if (pset == "gcsi") {
    anno$tissue <- tissue_meta$Tissue_supergroup[match(anno$Cell_line, tissue_meta$Cell_line)]
  } else if (pset == "ccle") {
    anno$tissue <- str_to_title(gsub("_", " ", tissue_meta$tissue[match(anno$Cell_Line, tissue_meta$Cell_Line)]))
  } else {
    anno$tissue <- str_to_title(tissue_meta$Factor.Value.organism.part.[match(anno$sample_alias, tissue_meta$Source.Name)])
  }
  return(anno)
}


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

#' Filter transcripts with 0 expression across intersected cell lines
#' Used after subsetting to 48 common cell lines
#' 
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