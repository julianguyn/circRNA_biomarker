#' helper function to parse and format dataframes
parse_df <- function(df, label, cells) {
    df <- df[df$dataset == label,]
    df$dataset <- NULL
    rownames(df) <- cells
    return(df)
}

#' Function to perform pairwise correlations
correlate_profiles <- function(gcsi_df, ccle_df, gdsc_df, circ = TRUE) {

    if (circ == TRUE) {
        # order cells
        gcsi_df <- gcsi_df[order(gcsi_df$sample), ]
        ccle_df <- ccle_df[order(ccle_df$sample), ]
        gdsc_df <- gdsc_df[order(gdsc_df$sample), ]

        # keep only common transcripts
        df <- bind_rows(
            gcsi = gcsi_df[, colnames(gcsi_df) %in% common_transcripts],
            ccle = ccle_df[, colnames(ccle_df) %in% common_transcripts],
            gdsc = gdsc_df[, colnames(gdsc_df) %in% common_transcripts],
            .id = "dataset"
            ) %>%
            mutate(across(-dataset, ~replace_na(as.double(.), 0)))

        # parse datasets
        gcsi_df <- parse_df(df, "gcsi", gcsi_df$sample)
        ccle_df <- parse_df(df, "ccle", ccle_df$sample)
        gdsc_df <- parse_df(df, "gdsc", gdsc_df$sample)
    }

    # initiate dataframe for storing results
    correlations <- data.frame(matrix(nrow=0, ncol=4))
    colnames(correlations) <- c("Cells", "gCSI & CCLE", "gCSI & GDSC2", "CCLE & GDSC2")

    # loop through each common transcript
    for (i in 1:nrow(gcsi_df)) {

        cell <- rownames(gcsi_df)[i]

        # compute correlations of transcript expression for pairs of psets
        gcsi_ccle <- suppressWarnings(cor(x = as.numeric(gcsi_df[i,]), y = as.numeric(ccle_df[i,]), method = "spearman")) #gCSI vs CCLE
        gcsi_gdsc <- suppressWarnings(cor(x = as.numeric(gcsi_df[i,]), y = as.numeric(gdsc_df[i,]), method = "spearman")) #gCSI vs GDSC
        ccle_gdsc <- suppressWarnings(cor(x = as.numeric(ccle_df[i,]), y = as.numeric(gdsc_df[i,]), method = "spearman")) #GDSC vs CCLE
    
        # combine results
        correlations <- rbind(correlations, 
                            data.frame(Cells = cell, gCSI_CCLE = gcsi_ccle, gCSI_GDSC = gcsi_gdsc, GDSC_CCLE = ccle_gdsc))
    }

    return(correlations)
}

#' Function to perform pairwise correlations
umap_cells <- function(gcsi_df, ccle_df, gdsc_df, circ = TRUE) {

    if (circ == TRUE) {
        # order cells
        gcsi_df <- gcsi_df[order(gcsi_df$sample), ]
        ccle_df <- ccle_df[order(ccle_df$sample), ]
        gdsc_df <- gdsc_df[order(gdsc_df$sample), ]

        # keep only common transcripts
        df <- bind_rows(
            gcsi = gcsi_df[, colnames(gcsi_df) %in% common_transcripts],
            ccle = ccle_df[, colnames(ccle_df) %in% common_transcripts],
            gdsc = gdsc_df[, colnames(gdsc_df) %in% common_transcripts],
            .id = "dataset"
            ) %>%
            mutate(across(-dataset, ~replace_na(as.double(.), 0)))

        # parse datasets
        gcsi_df <- parse_df(df, "gcsi", gcsi_df$sample)
        ccle_df <- parse_df(df, "ccle", ccle_df$sample)
        gdsc_df <- parse_df(df, "gdsc", gdsc_df$sample)
    }

    # save cell line labels
    cell_label <- c(gcsi_df$sample, ccle_df$sample, gdsc_df$sample)
    rownames(gcsi_df) <- NULL
    rownames(ccle_df) <- NULL
    rownames(gdsc_df) <- NULL

    # perform umap
    umap_df <- umap(rbind(gcsi_df, ccle_df, gdsc_df))
    umap_df <- umap_df$layout |> as.data.frame()
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$cell_line <- cell_label
    umap_df$dataset <- rep(names(pset_pal), each = 48)
    return(umap_df)
}