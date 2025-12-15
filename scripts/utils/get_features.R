#' Get transcript features from expression matrix
#' 
#' @param df data.frame. Expression matrix (transcripts as cols)
#' 
get_features <- function(df) {

    # initiate dataframe to store results
    res <- data.frame(matrix(nrow = ncol(df), ncol = 4))
    rownames(res) <- colnames(df)
    colnames(res) <- c("MedianExp", "Length", "NExons", "GC")

    # compute median expression
    res$MedianExp <- colSums(df)/48 |> as.numeric()

    # compute length
    coords <- str_split(rownames(res), "\\.", simplify = TRUE) |> as.data.frame()
    colnames(coords) <- c("chr", "start", "end")
    res$Length <- as.numeric(coords$end) - as.numeric(coords$start) + 1
 
    # compute GC context
    for (i in 1:nrow(coords)) {
        transcript_seq <- BSgenome::getSeq(
            genome, 
            names = coords$chr[i], 
            start = as.numeric(coords$start[i]), 
            end = as.numeric(coords$end[i])
        )
        gc_percent <- letterFrequency(transcript_seq, letters = c("G", "C"), as.prob = TRUE) |> sum() * 100
        res$GC[i] <- gc_percent
    }

    # create GRanges of transcript coords
    gr <- makeGRangesFromDataFrame(coords)
    mcols(gr)$transcriptID <- colnames(df)

    # extract number of exons
    for (i in seq_along(gr)) {
        transcript <- gr[i]
        overlaps <- findOverlaps(exons, transcript)
        hits <- exons[queryHits(overlaps)]
        unique_exons <- GRanges()

        # ensure at least one exon is contained in transcript
        for (gene in names(hits)) {
            indiv_exons <- hits[[gene]]
            olaps <- findOverlaps(indiv_exons, transcript, type = "within")
            if (length(olaps) > 0)  {
                unique_exons <- c(unique_exons, indiv_exons[queryHits(olaps)])
            }
        }
        unique_exons <- reduce(unique_exons)
        res$NExons[i] <- length(unique_exons)
    }
    return(res)
}

#' Format stability dataframes for Elastic net
#' 
#' @param stability dataframe. Output of si_Distribution
#' @param gcsi dataframe. Output of get_features()
#' @param ccle dataframe. Output of get_features()
#' @param gdsc dataframe. Output of get_features()
#' @param label string
#'
format_stability <- function(stability, gcsi, ccle, gdsc, label) {

    # reformat stability dataframe
    stability <- reshape2::dcast(stability, transcript ~ label, value.var = "stability")
    colnames(stability) <- c("transcript", "gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")
    
    # get median expression
    stability$gcsi_median <- stability$ccle_median <- stability$gdsc_median <- 0
    for (i in seq_along(stability$transcript)) {
        transcript <- stability$transcript[i]

        if (transcript %in% rownames(gcsi)) {
        stability$gcsi_median[i] <- gcsi$MedianExp[rownames(gcsi) == transcript]
        }
        if (transcript %in% rownames(ccle)) {
        stability$ccle_median[i] <- ccle$MedianExp[rownames(ccle) == transcript]
        }
        if (transcript %in% rownames(gdsc)) {
        stability$gdsc_median[i] <- gdsc$MedianExp[rownames(gdsc) == transcript]
        }
    }

    # get GC, nexons, and length
    feats <- data.frame(
        transcript = c(rownames(gcsi), rownames(ccle), rownames(gdsc)),
        gc = c(gcsi$GC, ccle$GC, gdsc$GC),
        nexons = c(gcsi$NExons, ccle$NExons, gdsc$NExons),
        length = c(gcsi$Length, ccle$Length, gdsc$Length)
    )
    feats <- unique(feats)
    feats <- feats[match(stability$transcript, feats$transcript),]
    
    stability$gc <- feats$gc
    stability$n_exon <- feats$nexons
    stability$length <- feats$length

    # save file
    filename <- paste0("../results/data/temp/", label, "_stability.csv")
    write.csv(stability, file = filename, quote = FALSE, row.names = FALSE)
} 