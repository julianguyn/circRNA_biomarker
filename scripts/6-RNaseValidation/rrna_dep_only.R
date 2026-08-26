# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

############################################################
# Load in circRNA expression data
############################################################
# from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113120

counts <- fread("data/GSE113120_rnaseq_counts_circRNA.txt", data.table = FALSE)
rownames(counts) <- counts$V1
counts$V1 <- NULL

############################################################
# Format counts matrix
############################################################

counts <- colSums(counts) |> as.data.frame()
colnames(counts) <- "Counts"
counts$sampleid <- gsub("_pos|_neg|plus|-", "", rownames(counts))
counts$RNaseR <- rep(c("RNase R-", "RNase R+"), 8)


ggplot(counts, aes(x = sampleid, y = Counts, fill = RNaseR)) +
    geom_bar(stat = "identity", position = "dodge")

p <- ggplot(counts, aes(x=sampleid, y=log2(Counts), fill=RNaseR)) + 
    geom_bar(stat="identity", width = 0.7, position = position_dodge(), color = "black") + 
    scale_y_continuous(limits = c(0, 18), expand=c(0,0)) +
    #scale_fill_manual("Method", values=protocol_pal, limits=c("polyA", "ribo0"), labels=c("poly(A)-selected", "rRNA-depleted")) + 
    theme_bw() +
    theme(strip.background = element_rect(fill = "white", color = "black"),
            legend.key.size = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), 
            axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45)) +
    labs(x = "Cell line \n", y = "Log2 Normalized Counts")
    