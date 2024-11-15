# Script to correlate circRNA expression across biological replicates (alternative to UMAP)

# load libraries
suppressMessages(library(data.table))

# load processed expression data
load("../results/data/umapdf.RData")

# function to parse processed dataframe and perform pairwise correlations
corr_reps <- function(df, cell_line_labels) {
    
    ## current df structure: 144 rows (48 cell lines x 3 PSets) x unique cell lines (columns)

    # split dataframe into three psets
    p1 <- df[1:48,]
    p2 <- df[49:96,]
    p3 <- df[97:144,]

    # add cell line labels
    rownames(p1) <- cell_line_labels[1:48]
    rownames(p2) <- cell_line_labels[49:96]
    rownames(p3) <- cell_line_labels[97:144]

    # order cell lines
    p1 <- p1[order(rownames(p1)),]
    p2 <- p2[order(rownames(p2)),]
    p3 <- p3[order(rownames(p3)),]

    # initiate dataframe for storing results
    correlations <- data.frame(matrix(nrow=0, ncol=4))
    colnames(correlations) <- c("Cells", "gCSI_CCLE", "gCSI_GDSC", "GDSC_CCLE")

    # loop through each common transcript
    for (i in 1:nrow(p1)) {

        cell <- rownames(p1)[i]

        # compute correlations of transcript expression for pairs of psets
        p1_p2 <- suppressWarnings(cor(x = as.numeric(p1[i,]), y = as.numeric(p2[i,]), method = "spearman")) #gCSI vs CCLE
        p1_p3 <- suppressWarnings(cor(x = as.numeric(p1[i,]), y = as.numeric(p3[i,]), method = "spearman")) #gCSI vs GDSC
        p2_p3 <- suppressWarnings(cor(x = as.numeric(p2[i,]), y = as.numeric(p3[i,]), method = "spearman")) #GDSC vs CCLE
    
        # combine results
        correlations <- rbind(correlations, 
                            data.frame(Cells = cell, gCSI_CCLE = p1_p2, gCSI_GDSC = p1_p3, GDSC_CCLE = p2_p3))
    }

    return(correlations)
}

# compute correlations 
gexpr_corr <- corr_reps(gexpr_df, cell_line_gexpr)
isoform_corr <- corr_reps(isoform_df, cell_line_isoforms)
ciri_corr <- corr_reps(ciri_df, cell_line_ciri)
circ_corr <- corr_reps(circ_df, cell_line_circ)
cfnd_corr <- corr_reps(cfnd_df, cell_line_cfnd)
fcrc_corr <- corr_reps(fcrc_df, cell_line_fcrc)

save(ciri_corr, circ_corr, cfnd_corr, fcrc_corr, gexpr_corr, isoform_corr, file = "../results/data/corr_expr.RData")
