#' Function to compute transcript stability indices for a given pair of datasets
#' 
#' @param df1 dataframe. Dataframe with samples as rownames
#' @param df2 dataframe. Dataframe with samples as rownames
#' @param label1 str. Label of df1
#' @param label2 str. Label of df2
#' @param dataset str. Dataset label
#' @param random boolean. TRUE for random sampling of sample names
#' @param iter int. Number of iterations (for random sampling)
#' 
compute_stability_index <- function(df1, df2, label1, label2, dataset, random = FALSE, iter = 1) {

    # make label
    label <- paste(label1, label2, sep = "/")
    print(paste("Computing stability index for", label))

    # common samples and transcripts
    samples <- intersect(rownames(df1), rownames(df2))
    transcripts <- intersect(colnames(df1), colnames(df2))
    print(paste("-----Number of common transcripts:", length(transcripts)))

    # catch case where there are no overlapping transcripts
    if (length(transcripts) == 0) {
        return(data.frame(matrix(nrow=0, ncol=4)))
    }

    # order dataframes
    df1 <- df1[match(samples, rownames(df1)), match(transcripts, colnames(df1))]
    df2 <- df2[match(samples, rownames(df2)), match(transcripts, colnames(df2))]

    # initialize vectors to store results
    n <- length(transcripts)
    t_names <- rep(transcripts, iter)
    si_res  <- numeric(n * iter)

    idx <- 1
    for (i in seq_len(iter)) {

        if (random) {
            # shuffle the rows (samples)
            df1 <- df1[sample(samples), ]
            df2 <- df2[sample(samples), ]
        }

        # vectorized spearman across all transcripts
        corr <- diag(cor(df1, df2, method = "spearman")) |> suppressWarnings()

        # store results
        si_res[idx:(idx + n - 1)] <- corr
        idx <- idx + n
    }

    correlations <- data.frame(
        transcript = t_names,
        stability  = si_res
    )
    correlations$label <- label
    correlations$dataset <- dataset
    correlations$random <- ifelse(random == TRUE, "Random", "NonRandom")

    return(correlations)
}