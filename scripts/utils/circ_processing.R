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