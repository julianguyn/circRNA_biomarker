# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(stringr)
    library(viridis)
    library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(123)

source("utils/palettes.R")

args <- "circ"
dir <- "../data/seqdepth_cellline/"

############################################################
# Load in circRNA expression data
############################################################

# helper function to load in data
load_counts <- function(dir, filename) {
    counts <- fread(paste0(dir, filename), data.table = FALSE)
    rownames(counts) <- counts$sample
    counts$sample <- NULL
    return(counts)
}

ciri_gcsi <- load_counts(dir, "CIRI2/ciri_gcsi_counts.tsv")
ciri_gdsc <- load_counts(dir, "CIRI2/ciri_gdsc_counts.tsv")
ciri_ccle <- load_counts(dir, "CIRI2/ciri_ccle_counts.tsv")

circ_gcsi <- load_counts(dir, "CIRCexplorer2/circ_gcsi_counts.tsv")
circ_gdsc <- load_counts(dir, "CIRCexplorer2/circ_gdsc_counts.tsv")
circ_ccle <- load_counts(dir, "CIRCexplorer2/circ_ccle_counts.tsv")

cfnd_gcsi <- load_counts(dir, "circRNA_finder/cfnd_gcsi_counts.tsv")
cfnd_gdsc <- load_counts(dir, "circRNA_finder/cfnd_gdsc_counts.tsv")
cfnd_ccle <- load_counts(dir, "circRNA_finder/cfnd_ccle_counts.tsv")

fcrc_gcsi <- load_counts(dir, "find_circ/fcrc_gcsi_counts.tsv")
fcrc_gdsc <- load_counts(dir, "find_circ/fcrc_gdsc_counts.tsv")
fcrc_ccle <- load_counts(dir, "find_circ/fcrc_ccle_counts.tsv")

############################################################
# Load in circRNA read counts
############################################################

gcsi_reads <- fread(paste0(dir, "readCounts/gcsi.tsv"), data.table = FALSE)
ccle_reads <- fread(paste0(dir, "readCounts/ccle.tsv"), data.table = FALSE)
gdsc_reads <- fread(paste0(dir, "readCounts/gdsc.tsv"), data.table = FALSE)

############################################################
# Count all transcripts
############################################################

# count reads
count_reads <- function(counts, reads, label) {
    count <- rowSums(counts) |> as.data.frame()
    colnames(count) <- "Count"
    count$Seq_Depth <- reads$avg_counts[match(rownames(count), reads$sample)]
    count$Label <- label
    count$Sample <- rownames(count)
    rownames(count) <- NULL
    return(count)
}

ciri_gcsi <- count_reads(ciri_gcsi, gcsi_reads, "ciri_gcsi")
ciri_gdsc <- count_reads(ciri_gdsc, gdsc_reads, "ciri_gdsc")
ciri_ccle <- count_reads(ciri_ccle, ccle_reads, "ciri_ccle")

circ_gcsi <- count_reads(circ_gcsi, gcsi_reads, "circ_gcsi")
circ_gdsc <- count_reads(circ_gdsc, gdsc_reads, "circ_gdsc")
circ_ccle <- count_reads(circ_ccle, ccle_reads, "circ_ccle")

cfnd_gcsi <- count_reads(cfnd_gcsi, gcsi_reads, "cfnd_gcsi")
cfnd_gdsc <- count_reads(cfnd_gdsc, gdsc_reads, "cfnd_gdsc")
cfnd_ccle <- count_reads(cfnd_ccle, ccle_reads, "cfnd_ccle")

fcrc_gcsi <- count_reads(fcrc_gcsi, gcsi_reads, "fcrc_gcsi")
fcrc_gdsc <- count_reads(fcrc_gdsc, gdsc_reads, "fcrc_gdsc")
fcrc_ccle <- count_reads(fcrc_ccle, ccle_reads, "fcrc_ccle")

save(
    ciri_gcsi, ciri_gdsc, ciri_ccle,
    circ_gcsi, circ_gdsc, circ_ccle,
    cfnd_gcsi, cfnd_gdsc, cfnd_ccle,
    fcrc_gcsi, fcrc_gdsc, fcrc_ccle,
    file = "../results/data/seqdepth_temp.RData"
)

############################################################
# Create toPlot
############################################################

toPlot <- rbind(
    ciri_gcsi, ciri_gdsc, ciri_ccle,
    circ_gcsi, circ_gdsc, circ_ccle,
    cfnd_gcsi, cfnd_gdsc, cfnd_ccle,
    fcrc_gcsi, fcrc_gdsc, fcrc_ccle
)
toPlot$Pipeline <- sub("_.*", "", toPlot$Label)
toPlot$PSet <- sub(".*_", "", toPlot$Label)

toPlot$Pipeline <- factor(toPlot$Pipeline, levels = c("ciri", "circ", "cfnd", "fcrc"), labels = names(pipeline_pal))
toPlot$PSet <- factor(toPlot$PSet, levels = c("gcsi", "ccle", "gdsc"), labels = names(pset_pal))



############################################################
# Plot
############################################################

plot_seq_depth <- function(toPlot, pset) {

    toPlot <- toPlot[toPlot$PSet == pset,]

    toPlot$Seq_Depth <- as.numeric(toPlot$Seq_Depth)
    toPlot <- toPlot[order(toPlot$Seq_Depth),]
    toPlot$Sample <- factor(toPlot$Sample, levels = unique(toPlot$Sample))

    p1 <- ggplot(toPlot, aes(x = Pipeline, y = Sample, fill = log2(Count))) +
        geom_tile() +
        scale_fill_viridis(limits=c(0, 16), option="mako", direction = -1) +
        theme_minimal() +
        ggtitle(pset)

    p2 <- ggplot(toPlot, aes(x = Seq_Depth, y = Sample)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(
            axis.title.y = element_blank()
        ) +
        xlim(c(0, 263170100)) +
        labs(x = "Sequencing Depth")

    p <- p1 + p2 + plot_layout(guides = "collect")
    return(p)
}

plot_seq_depth(toPlot, "gCSI")
plot_seq_depth(toPlot, "CCLE")
plot_seq_depth(toPlot, "GDSC2")

############################################################
# Create list of unique transcripts (for venn diagram)
############################################################

# create list object of transcripts
ciri_transcripts <- list(gCSI = colnames(ciri_gcsi), CCLE = colnames(ciri_ccle), GDSC2 = colnames(ciri_gdsc))
circ_transcripts <- list(gCSI = colnames(circ_gcsi), CCLE = colnames(circ_ccle), GDSC2 = colnames(circ_gdsc))
cfnd_transcripts <- list(gCSI = colnames(cfnd_gcsi), CCLE = colnames(cfnd_ccle), GDSC2 = colnames(cfnd_gdsc))
fcrc_transcripts <- list(gCSI = colnames(fcrc_gcsi), CCLE = colnames(fcrc_ccle), GDSC2 = colnames(fcrc_gdsc))

############################################################
# Create proportion plots
############################################################

toPlot <- rbind(
    count_prop(ciri_transcripts, "CIRI2"),
    count_prop(circ_transcripts, "CIRCexplorer2"),
    count_prop(cfnd_transcripts, "circRNA_finder"),
    count_prop(fcrc_transcripts, "find_circ")
)
toPlot$pipeline <- factor(toPlot$pipeline, levels = names(pipeline_pal))

# proportion bar plot
filename <- paste("../results/figures/figure1/proportion_pipelines_", analysis, ".png")
png(filename, width=160, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(fill = Var1, y = Freq, x = pipeline)) + 
  geom_bar(position = "fill", stat = "identity", color = "black") +
  #geom_text(aes(label = ifelse(Var1 %in% c("gCSI only", "CCLE only", "GDSC2 only"), label, "")), position = position_fill(vjust = 0.5)) +
  theme_classic() + 
  theme(legend.key.size = unit(0.5, 'cm')) +
  scale_fill_manual(values = prop_pal) +
  labs(fill = "Category", x = "Pipeline", y = "Proportion of Unique Transcripts")
dev.off()


############################################################
# Plot abundance per dataset / pipeline
############################################################

# minor formating
ciri_gdsc[is.na(ciri_gdsc)] <- 0
circ_gdsc[is.na(circ_gdsc)] <- 0
cfnd_gdsc[is.na(cfnd_gdsc)] <- 0
fcrc_gdsc[is.na(fcrc_gdsc)] <- 0

# create data frame of counts for plotting
df <- data.frame(
    Count = c(
        sum(ciri_gcsi), sum(ciri_ccle), sum(ciri_gdsc),
        sum(circ_gcsi), sum(circ_ccle), sum(circ_gdsc),
        sum(cfnd_gcsi), sum(cfnd_ccle), sum(cfnd_gdsc),
        sum(fcrc_gcsi), sum(fcrc_ccle), sum(fcrc_gdsc)
    ),
    PSet = c(rep(names(pset_pal), 4)),
    Pipeline = c(names(pipeline_pal, each = 3)))
df$Pipeline <- factor(df$Pipeline, levels = names(pipeline_pal))
df$PSet <- factor(df$PSet, levels = names(pset_pal))

# plot counts
filename <- paste0("../results/figures/figure1/counts_", analysis, ".png")
png(filename, width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(df, aes(x = Pipeline, y = log2(Count), fill = PSet)) +
    geom_bar(stat="identity", position = "dodge", color = "black") +
    scale_fill_manual(values = pset_pal) + 
    scale_y_continuous(limits = c(0, 13), expand=c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.4, 'cm')
    ) + 
    labs(y = "Log2 Normalized Counts")
dev.off()