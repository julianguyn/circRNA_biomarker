summarizeCIRI <- function(dir_path){
  
  circRNA_ids <- c() #holds all circRNA_ids for dataset
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T) 
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(ciri_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    #sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", ciri_files[f])
    sample_name <- gsub("*.tsv", "", sample_name)   
    sample <- sample[which(sample$circRNA_type=="exon"),] #keep only exonic circRNA's
    sample <- sample[which(sample$junction_reads >= 2),] #keep only circRNA's with >= 2 junction reads
    
    #convert coordinate from 1-based to 0-based to match CIRCexplorer2 coordinates
    sample$circRNA_start <- sample$circRNA_start - 1
    sample$circRNA_ID <- paste0(sample$chr, ":", sample$circRNA_start, "|", sample$circRNA_end)
    
    circRNA_ids <- c(sample$circRNA_ID, circRNA_ids)
    count <- sum(sample$junction_reads) #total number of junction reads per sample
    circ_counts_df[f,] <- c(sample_name,count)
  }
  #xx <- list(circ_counts_df, circRNA_ids)
  #names(xx)[1] <- "circ_counts_df"
  #names(xx)[2] <- "circRNA_ids"
  return(circ_counts_df)
} 

summarizeCIRCexplorer <- function(dir_path){
  
  circRNA_ids <- c() #holds all circRNA_ids for dataset
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = TRUE, recursive = TRUE)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(circ_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    
    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample <- sample[which(sample$circType=="circRNA"),] #keep only exonic circRNA's
    sample <- sample[which(sample$readNumber >= 2),] #keep only circRNA's with junction reads >= 2
    
    sample$circRNA_ID <- paste0(sample$chrom,":",sample$start,"|",sample$end) #create circRNA_id using chr, start & end coordinates
    circRNA_ids <- c(sample$circRNA_ID, circRNA_ids)
    count <- sum(sample$readNumber) #total # of junction reads per sample
    circ_counts_df[f,] <- c(sample_name,count)
  }
  #xx <- list(circ_counts_df, circRNA_ids)
  #names(xx)[1] <- "circ_counts_df"
  #names(xx)[2] <- "circRNA_ids"
  circ_counts_df$sample <- gsub("_circularRNA_known", "", circ_counts_df$sample)
  return(circ_counts_df)
} 


filterCIRI <- function(nonRNAseR_dir, RNAseR_dir, suff){
  #"*neg.tsv$" <- for ribominus
  #"\\.tsv$"<- for poly-A selected data
  ciri_files <- list.files(path=nonRNAseR_dir, pattern = suff, full.names = T)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(ciri_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    #sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", ciri_files[f])
    sample_name <- gsub("*.tsv", "", sample_name)
    match_name <- gsub("ccle|gcsi", "",ciri_files[f])
    #match_name <- gsub("\\..*","", match_name)
    match_name <- gsub(".*/","", match_name)
    match_name <- gsub("neg","", match_name)
    match_name <- gsub("*.tsv", "", match_name)
    
    hansen_match <- read.table(file = paste0(RNAseR_dir, match_name,"pos.tsv"), 
                               sep = '\t', 
                               skip = 1,
                               header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    colnames(hansen_match) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    
    sample_match <- sample[which(sample$circRNA_type=="exon"),] #keep exonic circRNA's only
    sample_match <- sample_match[which(sample_match$junction_reads >= 2),] #keep only circRNA's with junction reads >= 2
    
    hansen_match <- hansen_match[which(hansen_match$circRNA_type=="exon"),] #keep exonic circRNA's only
    hansen_match <- hansen_match[which(hansen_match$junction_reads >= 2),]  #keep only circRNA's with junction reads >= 2 
    
    sample_match <- sample_match[which(sample_match$circRNA_ID %in% hansen_match$circRNA_ID),]
    hansen_match <- hansen_match[which(hansen_match$circRNA_ID %in% sample_match$circRNA_ID),]
    
    filtered <- sample_match[hansen_match$junction_reads/sample_match$junction_reads >= 1.5,] #keep only circRNA's with at least 1.5-fold enrichment in matching RNAse-R sample
    
    count <- sum(filtered$junction_reads) #total # of junction reads per sample
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 


filterCIRC <- function(nonRNAseR_dir, RNAseR_dir, suff){
  #"*neg*" <- for ribominus
  #"\\.txt$"<- for poly-A selected data
  circ_files <- list.files(path=nonRNAseR_dir,  pattern = suff, full.names = T, recursive = TRUE)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(circ_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    
    sample_name <- gsub("_circularRNA_known.txt","", circ_files[f])
    sample_name <- gsub(".*/","", sample_name)
    match_name <- gsub("neg", "",sample_name)
    match_name <- gsub("ccle|gcsi", "",match_name)
    
    hansen_match <- read.table(file = paste0(RNAseR_dir, "/" ,match_name, "pos", "/", match_name,"pos", "_circularRNA_known.txt"),
                               sep = '\t',
                               header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    colnames(hansen_match) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    
    sample_match <- sample[which(sample$circType=="circRNA"),]
    sample_match <- sample_match[which(sample_match$readNumber >= 2),] #keep only circRNA's with junction reads >= 2
    
    hansen_match <- hansen_match[which(hansen_match$circType=="circRNA"),]
    hansen_match <- hansen_match[which(hansen_match$readNumber >= 2),]
    
    sample_match$pos <- paste0(sample_match$chrom,sample_match$start,sample_match$end)
    hansen_match$pos <- paste0(hansen_match$chrom,hansen_match$start,hansen_match$end)
    
    sample_match <- sample_match[which(sample_match$pos %in% hansen_match$pos),]
    hansen_match <- hansen_match[which(hansen_match$pos %in% sample_match$pos),]
    hansen_match  <- hansen_match[match(sample_match$pos, hansen_match$pos),]
    
    filtered <- sample_match[hansen_match$readNumber/sample_match$readNumber >= 1.5,] #keep only circRNA's with at least 1.5-fold enrichment in matching RNAse-R sample
    
    count <- sum(filtered$readNumber)
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 
