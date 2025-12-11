#' Function to count proportions of unique transcript
#' across overlaps of PSets
#' 
#' @param transcripts list. List object of transcripts
#' @param pipeline str. Detection pipeline label
#' 
count_prop <- function(transcripts, pipeline) {
    
    gcsi <- transcripts$gCSI
    ccle <- transcripts$CCLE
    gdsc <- transcripts$GDSC2

    # get all transcripts
    all_transcripts <- unique(c(gcsi, ccle, gdsc))

    # label transcript detection per dataset
    transcript_df <- data.frame(
        transcript = all_transcripts,
        gcsi = all_transcripts %in% gcsi,
        ccle = all_transcripts %in% ccle,
        gdsc = all_transcripts %in% gdsc
    )

    # label transcript detection
    transcript_df$category <- with(transcript_df, 
        ifelse(gcsi & !ccle & !gdsc, "gCSI only",
        ifelse(!gcsi & ccle & !gdsc, "CCLE only",
        ifelse(!gcsi & !ccle & gdsc, "GDSC2 only",
        ifelse(gcsi & ccle & !gdsc, "gCSI & CCLE",
        ifelse(gcsi & !ccle & gdsc, "gCSI & GDSC2",
        ifelse(!gcsi & ccle & gdsc, "CCLE & GDSC2",
        "All datasets"))))))
    )

    # create dataframe for plotting
    toPlot <- table(transcript_df$category) |> as.data.frame()
    toPlot$Var1 <- factor(
        toPlot$Var1, 
        levels = c(
            "All datasets",
            "gCSI & CCLE",
            "gCSI & GDSC2",
            "CCLE & GDSC2",
            "gCSI only",
            "CCLE only",
            "GDSC2 only"
        ))
    toPlot$Prop <- round(toPlot$Freq / sum(toPlot$Freq) * 100, digits = 2)
    toPlot$label <- paste0(toPlot$Freq, "\n(", toPlot$Prop, "%)")
    toPlot$pipeline <- pipeline

    return(toPlot)
}