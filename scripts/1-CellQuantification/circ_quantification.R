# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(ggvenn)
    library(stringr)
    library(ComplexHeatmap)
})

options(stringsAsFactors = FALSE)
set.seed(123)

source("utils/palettes.R")
source("utils/transcript_quantification.R")

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("circ", "GE")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

############################################################
# Load in circRNA expression data
############################################################

# get input dirs
dir <- switch(
  analysis,
  circ = "../data/processed_cellline/common_samples/",
  GE = "../data/processed_cellline/GE_common_samples/"
)

# helper function to load in data
load_counts <- function(dir, filename) {
    counts <- fread(paste0(dir, filename), data.table = FALSE)
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
# Create list of unique transcripts (for venn diagram)
############################################################

# create list object of transcripts
ciri_transcripts <- list(gCSI = colnames(ciri_gcsi), CCLE = colnames(ciri_ccle), GDSC2 = colnames(ciri_gdsc))
circ_transcripts <- list(gCSI = colnames(circ_gcsi), CCLE = colnames(circ_ccle), GDSC2 = colnames(circ_gdsc))
cfnd_transcripts <- list(gCSI = colnames(cfnd_gcsi), CCLE = colnames(cfnd_ccle), GDSC2 = colnames(cfnd_gdsc))
fcrc_transcripts <- list(gCSI = colnames(fcrc_gcsi), CCLE = colnames(fcrc_ccle), GDSC2 = colnames(fcrc_gdsc))

# create list object of transcripts across all pipelines
all_comparison <- list(
    CIRI2 = c(colnames(ciri_gcsi), colnames(ciri_ccle), colnames(ciri_gdsc)), 
    CIRCexplorer2 = c(colnames(circ_gcsi), colnames(circ_ccle), colnames(circ_gdsc)), 
    circRNA_finder = c(colnames(cfnd_gcsi), colnames(cfnd_ccle), colnames(cfnd_gdsc)),
    find_circ = c(colnames(fcrc_gcsi), colnames(fcrc_ccle), colnames(fcrc_gdsc))
)

############################################################
# Plot venn diagram
############################################################

# helper function to plot venn diagram
plot_venn <- function(list_obj, pal, title = "") {
    p <- ggvenn(list_obj, fill_color = pal, stroke_size = 0.5, set_name_size = 4) +
        theme(plot.title = element_text(hjust = 0.5, size = 15)) +
        labs(title = title)
    return(p)
}

# plot transcript detection per pipeline
p1 <- plot_venn(ciri_transcripts, pset_pal, "CIRI2")
p2 <- plot_venn(circ_transcripts, pset_pal, "CIRCexplorer2")
p3 <- plot_venn(cfnd_transcripts, pset_pal, "circRNA_finder")
p4 <- plot_venn(fcrc_transcripts, pset_pal, "find_circ")

filename <- paste0("../rsults/figures/suppfigure1/indiv_venn_", analysis, ",png")
png(filename, width=500, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = FALSE)
dev.off()

# plot transcript detection across all pipelines
filename <- paste0("../rsults/figures/suppfigure1/pipeline_venn_", analysis, ",png")
png(filename, width=200, height=150, units='mm', res = 600, pointsize=80)
plot_venn(all_comparison, pipeline_pal)
dev.off()


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
# Create Upset plot
############################################################

# create list object of transcripts for upset plot
toPlot <- make_comb_mat(list(
    gCSI_CIRI2 = colnames(ciri_gcsi),
    CCLE_CIRI2 = colnames(ciri_ccle),
    GDSC_CIRI2 = colnames(ciri_gdsc),
    gCSI_CIRCexplorer2 = colnames(circ_gcsi),
    CCLE_CIRCexplorer2 = colnames(circ_ccle),
    GDSC_CIRCexplorer2 = colnames(circ_gdsc),
    gCSI_circRNA_finder = colnames(cfnd_gcsi),
    CCLE_circRNA_finder = colnames(cfnd_ccle),
    GDSC_circRNA_finder = colnames(cfnd_gdsc),
    gCSI_find_circ = colnames(fcrc_gcsi),
    CCLE_find_circ = colnames(fcrc_ccle),
    GDSC_find_circ = colnames(fcrc_gdsc)
))

# remove combinations of less than 100 pairs
toPlot <- toPlot[comb_size(toPlot) >= 100]

# specify set orders
set_order <- c(
    "gCSI_CIRI2",
    "CCLE_CIRI2",
    "GDSC_CIRI2",
    "gCSI_CIRCexplorer2",
    "CCLE_CIRCexplorer2",
    "GDSC_CIRCexplorer2",
    "gCSI_circRNA_finder",
    "CCLE_circRNA_finder",
    "GDSC_circRNA_finder",
    "gCSI_find_circ",
    "CCLE_find_circ",
    "GDSC_find_circ"
)

# upset plot
filename <- paste0("../results/figures/suppfigure1/upset_", anlaysis, ".pdf")
pdf(filename, width=10, height=5)
UpSet(toPlot, 
    set_order = set_order,
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
    comb_order = order(comb_size(toPlot)))
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