# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(viridis)
    library(reshape2)
})

options(stringsAsFactors = F)


############################################################
# Load in data 
############################################################

# analysis = c("circ", "GE")  

if (analysis == "circ") {
    load("../data/processed_lung/circ_lung_expression.RData") # from circ_lung.R
} else if (analysis == "GE") {
    ciri_polyA <- fread("../data/processed_lung/GE/ciri_polyA_counts.tsv")
    ciri_ribo0 <- fread("../data/processed_lung/GE/ciri_ribo0_counts.tsv")
    circ_polyA <- fread("../data/processed_lung/GE/circ_polyA_counts.tsv")
    circ_ribo0 <- fread("../data/processed_lung/GE/circ_ribo0_counts.tsv")
    cfnd_polyA <- fread("../data/processed_lung/GE/cfnd_polyA_counts.tsv")
    cfnd_ribo0 <- fread("../data/processed_lung/GE/cfnd_ribo0_counts.tsv")
    fcrc_polyA <- fread("../data/processed_lung/GE/fcrc_polyA_counts.tsv")
    fcrc_ribo0 <- fread("../data/processed_lung/GE/fcrc_ribo0_counts.tsv")
}

# use the lung metadata to find the polyA samples to keep
polyA_meta <- read.csv("../data/rnaseq_meta/lung_polyA.tsv", sep = "\t")
ribo0_meta <- read.csv("../data/rnaseq_meta/lung_ribozero.csv")

# get library sizes
#ribo0_lib <- read.table("../data/processed_circRNA/ribo0.out")
#ribo0_anno <- read.table("../data/processed_circRNA/ribo0.txt")
#write.table(ribo0_anno, file = "../data/processed_circRNA/rib0_lib.tsv", quote = F, sep = "\t", col.names = T)
polyA_lib <- read.table("../data/processed_circRNA/polyA_lib.tsv")
ribo0_lib <- read.table("../data/processed_circRNA/ribo0_lib.tsv")

# function to average read counts for paired-end samples
avg_reads <- function(lib) {

  for (sample in unique(lib$sample)) {
    idx <- grep(sample, lib$sample)
    #avg read counts
    avg <- mean(lib$libsize[idx])
    #rm paired-end values and replace with average
    lib <- lib[-idx,]
    lib <- rbind(lib, data.frame(sample = sample, libsize = avg))
  }

  return(lib)
}
polyA_lib <- avg_reads(polyA_lib)
polyA_lib$match <- ribo0_meta[match(polyA_meta[match(polyA_lib$sample, polyA_meta$sample),]$TB_id, ribo0_meta$TB_id),]$helab_id
ribo0_lib <- avg_reads(ribo0_lib)

# format polyA sample names
ciri_polyA <- ciri_polyA[rownames(ciri_polyA) %in% polyA_meta$sample,]
rownames(ciri_polyA) <- ribo0_meta[match(polyA_meta[match(rownames(ciri_polyA), polyA_meta$sample),]$TB_id, ribo0_meta$TB_id),]$helab_id
circ_polyA <- circ_polyA[gsub("_circularRNA_known", "", rownames(circ_polyA)) %in% polyA_meta$sample,]
rownames(circ_polyA) <- ribo0_meta[match(polyA_meta[match(gsub("_circularRNA_known", "", rownames(circ_polyA)), polyA_meta$sample),]$TB_id, ribo0_meta$TB_id),]$helab_id
cfnd_polyA <- cfnd_polyA[gsub("polyA/", "", gsub("filteredJunctions.bed", "", rownames(cfnd_polyA))) %in% polyA_meta$sample,]
rownames(cfnd_polyA) <- ribo0_meta[match(polyA_meta[match(gsub("polyA/", "", gsub("filteredJunctions.bed", "", rownames(cfnd_polyA))), polyA_meta$sample),]$TB_id, ribo0_meta$TB_id),]$helab_id
fcrc_polyA <- fcrc_polyA[gsub("polyA/", "", gsub("_circ.bed", "", rownames(fcrc_polyA))) %in% polyA_meta$sample,]
rownames(fcrc_polyA) <- ribo0_meta[match(polyA_meta[match(gsub("polyA/", "", gsub("_circ.bed", "", rownames(fcrc_polyA))), polyA_meta$sample),]$TB_id, ribo0_meta$TB_id),]$helab_id

# format ribo0 sample names
rownames(ciri_ribo0) <- gsub("_INPUT", "", rownames(ciri_ribo0))
rownames(circ_ribo0) <- gsub("_INPUT_circularRNA_known", "", rownames(circ_ribo0))
rownames(cfnd_ribo0) <- gsub("_INPUTfilteredJunctions.bed", "", gsub("ribo0/", "", rownames(cfnd_ribo0)))
rownames(fcrc_ribo0) <- gsub("_INPUT_circ.bed", "", gsub("ribo0/", "", rownames(ciri_ribo0)))


# function to normalize counts
cpm <- function(df, lib, seq) {

  # keep only chromosomal circRNAs
  rm <- colnames(df)[-grep("chr", colnames(df))]
  if(length(rm) > 0) {df <- df[,-which(colnames(df) %in% rm)]}

  # order the samples and match order with library size
  df <- df[order(rownames(df)),]
  ifelse(seq == "polyA", lib <- lib[match(rownames(df), lib$match),], lib <- lib[match(rownames(df), lib$sample),])

  # divide counts in each row by the corresponding library size for that sample
  df_cpm <- as.data.frame(lapply(df, function(col) col / lib$libsize * 1e6))

  # set rownames
  rownames(df_cpm) <- paste0("tumor",1:51)

  return(df_cpm)
}

ciri_polyA <- cpm(ciri_polyA, polyA_lib, "polyA")
circ_polyA <- cpm(circ_polyA, polyA_lib, "polyA")
cfnd_polyA <- cpm(cfnd_polyA, polyA_lib, "polyA")
fcrc_polyA <- cpm(fcrc_polyA, polyA_lib, "polyA")
ciri_ribo0 <- cpm(ciri_ribo0, ribo0_lib, "ribo0")
circ_ribo0 <- cpm(circ_ribo0, ribo0_lib, "ribo0")
cfnd_ribo0 <- cpm(cfnd_ribo0, ribo0_lib, "ribo0")
fcrc_ribo0 <- cpm(fcrc_ribo0, ribo0_lib, "ribo0")

save(ciri_polyA, circ_polyA, cfnd_polyA, fcrc_polyA, ciri_ribo0, circ_ribo0, cfnd_ribo0, fcrc_ribo0, file = "../results/data/circ_lung_expression.RData")


# merge data per pipeline for plotting
ciri_combined <- data.frame(samples = paste0("tumor",1:51), count_polyA = rowSums(ciri_polyA), count_ribo0 = rowSums(ciri_ribo0))
circ_combined <- data.frame(samples = paste0("tumor",1:51), count_polyA = rowSums(circ_polyA), count_ribo0 = rowSums(circ_ribo0))
cfnd_combined <- data.frame(samples = paste0("tumor",1:51), count_polyA = rowSums(cfnd_polyA), count_ribo0 = rowSums(cfnd_ribo0))
fcrc_combined <- data.frame(samples = paste0("tumor",1:51), count_polyA = rowSums(fcrc_polyA), count_ribo0 = rowSums(fcrc_ribo0))

# add labels
tumourID <- c(rownames(ciri_combined), rownames(circ_combined), rownames(cfnd_combined), rownames(fcrc_combined))
method <- c(rep("CIRI2", 51), rep("CIRCexplorer2", 51), rep("circRNA_finder", 51), rep("find_circ", 51))
lung_combined <- rbind(ciri_combined, circ_combined, cfnd_combined, fcrc_combined)
lung_combined$tumourID <- tumourID
lung_combined$method <- method
lung_combined <- melt(lung_combined)
lung_combined$variable <- factor(lung_combined$variable, levels = c("count_ribo0", "count_polyA"))
lung_combined$value <- log2(lung_combined$value + 1)

# plot heatmap of counts
png("../results/figures/C_heatmap.png", width=250, height=200, units='mm', res = 600, pointsize=80)
ggplot(lung_combined, aes(x = variable, y = tumourID, fill = value)) + 
  geom_tile() + facet_grid(~factor(method, levels=c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))) +
  theme_void() +
  theme(text = element_text(size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17, angle = 90),
        strip.text.x = element_text(size = 17),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  scale_fill_viridis(limits=c(0, 13), option="mako", direction = -1) +
  scale_x_discrete(labels = c("rRNA-\ndepleted", "poly(A)\nselected")) +
  guides(fill = guide_colourbar(barwidth = 0.2, barheight = 2, title = "log2 Normalized \ncircRNA Counts")) +
  labs(x = "\nMethod", y = "Tumour Sample")
dev.off()


# function to compute counts for each circRNA transcript per pipeline ############this needs to be sent to H4H
sum_transcripts <- function(polyA_df, ribo0_df) {
  
  # format data frames 
  polyA_df <- data.frame(t(polyA_df))
  ribo0_df <- data.frame(t(ribo0_df))
  
  polyA_df$count <- rowSums(polyA_df)
  ribo0_df$count <- rowSums(ribo0_df)

  polyA_df$transcript <- rownames(polyA_df)
  ribo0_df$transcript <- rownames(ribo0_df)
  
  # get common circRNA transcript IDs
  transcripts <- intersect(rownames(polyA_df), rownames(ribo0_df))
  polyA_df <- polyA_df[rownames(polyA_df) %in% transcripts,which(colnames(polyA_df) %in% c("count", "transcript"))]
  ribo0_df <- ribo0_df[rownames(ribo0_df) %in% transcripts,which(colnames(ribo0_df) %in% c("count", "transcript"))]

  # match order of transcripts
  polyA_df <- polyA_df[match(transcripts, polyA_df$transcript),]
  ribo0_df <- ribo0_df[match(transcripts, ribo0_df$transcript),]

  # create dataframe
  df <- data.frame(transcript = transcripts, polyA_count = polyA_df$count, ribo0_count = ribo0_df$count)

  # format data frame for subsequent plotting
  df <- melt(df)
  colnames(df) <- c("Transcript", "Method", "Count")

  # log2 + 1 normalize circRNA counts 
  df$Count <- log2(df$Count + 1)
  
  return(df)
}

#compute counts
ciri_counts <- sum_transcripts(ciri_polyA, ciri_ribo0)
circ_counts <- sum_transcripts(circ_polyA, circ_ribo0)
cfnd_counts <- sum_transcripts(cfnd_polyA, cfnd_ribo0)
fcrc_counts <- sum_transcripts(fcrc_polyA, fcrc_ribo0)

#function to plot density plots of circRNA counts
plot_density <- function(counts_df, title) {
  p <- ggplot(data = count_df) +
          geom_density(aes(x = Count, fill = Method), color = "black", alpha = 0.5) +
          scale_fill_manual(values=c("#4CC5AB", "#392C57")) +
          theme_classic() +
          theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), text = element_text(size = 15), 
                legend.key.size = unit(1, 'cm'), legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 17), 
                axis.text.x = element_text(size=12, vjust = 0.5), axis.text.y = element_text(size=15)) +
          labs(title = "CIRI2 Transcript Counts", y = "Density", x = "log2 Normalized Counts") +
          guides(fill = guide_legend(title = "Method", ncol = 2))
  return(p)
}

p1 <- plot_density(CIRI_counts, "CIRI2 Transcript Counts")
p2 <- plot_density(CIRC_counts, "CIRCexplorer2 Transcript Counts")
p3 <- plot_density(CFND_counts, "circRNA_finder Transcript Counts")
p4 <- plot_density(FCRC_counts, "find_circ Transcript Counts")

png("results/figures/B_density.png", width=275, height=400, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 1, common.legend = TRUE, legend = "bottom")
dev.off()



# function to compute proportion of unique transcripts per pipeline
plot_unique_transcripts <- function(polyA_df, ribo0_df, title) {

  # get transcripts detected from both seq protocols
  intersected_transcripts <- intersect(colnames(polyA_df), colnames(ribo0_df))

  # record number of unique transcripts and number of transcripts detected by both polyA and ribozero
  unique_transcripts <- data.frame(category = c("poly(A)", "RiboZero", "Both"),
                            count = c(length(colnames(polyA_df)) - length(intersected_transcripts),
                                      length(colnames(ribo0_df)) - length(intersected_transcripts),
                                      length(intersected_transcripts)))

  # compute percentages
  unique_transcripts$fraction = unique_transcripts$count / sum(unique_transcripts$count)

  # compute the cumulative percentages (top of each rectangle)
  unique_transcripts$ymax <- cumsum(unique_transcripts$fraction)
  unique_transcripts$ymin <- c(0, head(unique_transcripts$ymax, n = -1))

  p <- ggplot(unique_transcripts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
          geom_rect(color = "black") + coord_polar(theta="y") + 
          guides(fill = guide_legend(title = "Proportion of \nDetected Transcripts", ncol = 1)) + #make legend one column instead of 3
          scale_fill_manual(labels=c("Both poly(A)-selected and rRNA-depleted", "poly(A)-selected only", "rRNA-depleted only"), values = c("#3670A0", "#4CC5AB", "#392C57")) + 
          xlim(c(2, 4)) + theme_classic() +
          theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), text = element_text(size = 15), 
                legend.key.size = unit(1, 'cm'), legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 17), 
                axis.text.x = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          labs(title = title)
  
  return(p)
}

p1 <- plot_unique_transcripts(ciri_polyA, ciri_ribo0, "CIRI2")
p2 <- plot_unique_transcripts(circ_polyA, circ_ribo0, "CIRCexplorer2")
p3 <- plot_unique_transcripts(cfnd_polyA, cfnd_ribo0, "circRNA_finder")
p4 <- plot_unique_transcripts(fcrc_polyA, fcrc_ribo0, "find_circ")


png("../results/figures/C_proportion_unique.png", width=350, height=125, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = "bottom")
dev.off()


#### look at concordance of overlap between polyA and ribo0
plot_unique_concordance <- function(polyA_df, ribo0_df, label) {

  # get transcripts detected from both seq protocols
  intersected_transcripts <- intersect(colnames(polyA_df), colnames(ribo0_df))

  polyA_df <- polyA_df[,colnames(polyA_df) %in% intersected_transcripts]
  ribo0_df <- ribo0_df[,colnames(ribo0_df) %in% intersected_transcripts]

  # formating for comparison
  polyA_df$tumor <- rownames(polyA_df)
  ribo0_df$tumor <- rownames(ribo0_df)

  polyA_df <- melt(polyA_df)
  ribo0_df <- melt(ribo0_df)

  polyA_df$pairs <- paste0(polyA_df$tumor, polyA_df$variable)
  ribo0_df$pairs <- paste0(ribo0_df$tumor, ribo0_df$variable)

  # create merged dataframe
  toPlot <- data.frame(pair = polyA_df$pairs, polyA = polyA_df$value,
                       ribo0 = ribo0_df[match(polyA_df$pairs, ribo0_df$pairs),]$value)
  toPlot <- toPlot[-which(toPlot$ribo0 == 0 | toPlot$polyA == 0),]

  # compute difference between ribo0 and polyA
  toPlot$diff <- toPlot$ribo0 - toPlot$polyA
  toPlot <- toPlot[order(toPlot$diff),]
  toPlot$rank <- 1:nrow(toPlot)

  #waterfall plot
  #ggplot(toPlot, aes(x = rank, y = diff)) + geom_col(aes(fill = ifelse(diff > 0, "ribo0>polyA", "polyA>ribo0")), color = "black") +  theme_classic() + theme(axis.text.x = element_blank()) + labs(fill = "", x = "Tumor-circRNA pair", y = "Difference")

  greater <- nrow(toPlot[toPlot$diff > 5,])
  concord <- nrow(toPlot[toPlot$diff < 5 & toPlot$diff > -5,])
  lessthn <- nrow(toPlot[toPlot$diff < -5,])

  toPlot <- data.frame(difference = c("ribo0>polyA", "concordant", "polyA>ribo0"),
                       value = c(greater, concord, lessthn), label = label)

  return(toPlot)
}

df <- rbind(plot_unique_concordance(ciri_polyA, ciri_ribo0, "CIRI2"),
            plot_unique_concordance(circ_polyA, circ_ribo0, "CIRCexplorer2"), 
            plot_unique_concordance(cfnd_polyA, cfnd_ribo0, "circRNA_finder"),
            plot_unique_concordance(fcrc_polyA, fcrc_ribo0, "find_circ"))

df$label <- factor(df$label, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
png("PipelineComp.png", width = 8, height = 6, res = 600, units = "in")
ggplot(df, aes(x = label, y = value, fill = difference)) + geom_bar(position="fill", stat="identity") +  
      theme_classic() +
      scale_fill_manual(values = c("#646881", "#62BEC1", "#63595C")) +
      labs(fill = "Difference", y = "Porportion", x = "Pipeline")
dev.off()