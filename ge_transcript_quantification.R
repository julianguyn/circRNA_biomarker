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


############################################################
# Load in data 
############################################################

# load circRNA expression data
path <- "../data/processed_cellline/GE_common_samples/"

ciri_gcsi_ge <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F) 
ciri_gdsc_ge <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F) 
ciri_ccle_ge <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F) 

circ_gcsi_ge <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F) 
circ_gdsc_ge <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F) 
circ_ccle_ge <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F) 

cfnd_gcsi_ge <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F) 
cfnd_gdsc_ge <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F) 
cfnd_ccle_ge <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F) 

fcrc_gcsi_ge <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F) 
fcrc_gdsc_ge <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F) 
fcrc_ccle_ge <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F) 


############################################################
# Format exp dataframes for plotting
############################################################

# set palette for plotting
pal1 = c("#51C7AD", "#392C57", "#3670A0")
pal2 = c("#E8F6B1", "#A5DBB7", "#2088BC", "#26479D")

# remove sample names column for downstream analysis
ciri_gcsi_ge <- ciri_gcsi_ge[,-which(colnames(ciri_gcsi_ge) %in% c("sample"))]
ciri_ccle_ge <- ciri_ccle_ge[,-which(colnames(ciri_ccle_ge) %in% c("sample"))]
ciri_gdsc_ge <- ciri_gdsc_ge[,-which(colnames(ciri_gdsc_ge) %in% c("sample"))]

circ_gcsi_ge <- circ_gcsi_ge[,-which(colnames(circ_gcsi_ge) %in% c("sample"))]
circ_ccle_ge <- circ_ccle_ge[,-which(colnames(circ_ccle_ge) %in% c("sample"))]
circ_gdsc_ge <- circ_gdsc_ge[,-which(colnames(circ_gdsc_ge) %in% c("sample"))]

cfnd_gcsi_ge <- cfnd_gcsi_ge[,-which(colnames(cfnd_gcsi_ge) %in% c("sample"))]
cfnd_ccle_ge <- cfnd_ccle_ge[,-which(colnames(cfnd_ccle_ge) %in% c("sample"))]
cfnd_gdsc_ge <- cfnd_gdsc_ge[,-which(colnames(cfnd_gdsc_ge) %in% c("sample"))]

fcrc_gcsi_ge <- fcrc_gcsi_ge[,-which(colnames(fcrc_gcsi_ge) %in% c("sample"))]
fcrc_ccle_ge <- fcrc_ccle_ge[,-which(colnames(fcrc_ccle_ge) %in% c("sample"))]
fcrc_gdsc_ge <- fcrc_gdsc_ge[,-which(colnames(fcrc_gdsc_ge) %in% c("sample"))]



############################################################
# Figure 4: Venn diagrams of unique transcript detection
############################################################

# function to plot venn diagram
plot_venn <- function(list_obj, pal, title = "") {
    p <- ggvenn(list_obj, fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = title)
    return(p)
}

# create list object of transcripts
ciri_transcripts <- list(gCSI = colnames(ciri_gcsi_ge), CCLE = colnames(ciri_ccle_ge), GDSC2 = colnames(ciri_gdsc_ge))
circ_transcripts <- list(gCSI = colnames(circ_gcsi_ge), CCLE = colnames(circ_ccle_ge), GDSC2 = colnames(circ_gdsc_ge))
cfnd_transcripts <- list(gCSI = colnames(cfnd_gcsi_ge), CCLE = colnames(cfnd_ccle_ge), GDSC2 = colnames(cfnd_gdsc_ge))
fcrc_transcripts <- list(gCSI = colnames(fcrc_gcsi_ge), CCLE = colnames(fcrc_ccle_ge), GDSC2 = colnames(fcrc_gdsc_ge))

# create list object of transcripts across all pipelines
all_comparison <- list(CIRI2 = c(colnames(ciri_gcsi_ge), colnames(ciri_ccle_ge), colnames(ciri_gdsc_ge)), 
                       CIRCexplorer2 = c(colnames(circ_gcsi_ge), colnames(circ_ccle_ge), colnames(circ_gdsc_ge)), 
                       circRNA_finder = c(colnames(cfnd_gcsi_ge), colnames(cfnd_ccle_ge), colnames(cfnd_gdsc_ge)),
                       find_circ = c(colnames(fcrc_gcsi_ge), colnames(fcrc_ccle_ge), colnames(fcrc_gdsc_ge)))


# plot transcript detection per pipeline
p1 <- plot_venn(ciri_transcripts, pal1, "CIRI2")
p2 <- plot_venn(circ_transcripts, pal1, "CIRCexplorer2")
p3 <- plot_venn(cfnd_transcripts, pal1, "circRNA_finder")
p4 <- plot_venn(fcrc_transcripts, pal1, "find_circ")

png("../results/figures/figure4/venndiagram_per_pipeline_ge.png", width=500, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = FALSE)
dev.off()

# plot transcript detection across all pipelines
png("../results/figures/figure4/venndiagram_ge.png", width=200, height=150, units='mm', res = 600, pointsize=80)
plot_venn(all_comparison, pal2)
dev.off()


############################################################
# Figure 4: Upset plot of overlapping transcripts
############################################################

# function to create upset plot
plot_upset <- function(comb_mat, set_order) {
    p <- UpSet(comb_mat, set_order = set_order,
        top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
        comb_order = order(comb_size(comb_mat)))
    return(p)
}

# create list object of transcripts for upset plot
toPlot <- make_comb_mat(list(
            gCSI_CIRI2 = colnames(ciri_gcsi_ge), CCLE_CIRI2 = colnames(ciri_ccle_ge), GDSC_CIRI2 = colnames(ciri_gdsc_ge),
            gCSI_CIRCexplorer2 = colnames(circ_gcsi_ge), CCLE_CIRCexplorer2 = colnames(circ_ccle_ge), GDSC_CIRCexplorer2 = colnames(circ_gdsc_ge),
            gCSI_circRNA_finder = colnames(cfnd_gcsi_ge), CCLE_circRNA_finder = colnames(cfnd_ccle_ge), GDSC_circRNA_finder = colnames(cfnd_gdsc_ge),
            gCSI_find_circ = colnames(fcrc_gcsi_ge), CCLE_find_circ = colnames(fcrc_ccle_ge), GDSC_find_circ = colnames(fcrc_gdsc_ge)))

# remove combinations of less than 100 pairs
toPlot <- toPlot[comb_size(toPlot) >= 100]

# specify set orders
set_order1 = c("gCSI_CIRI2", "CCLE_CIRI2", "GDSC_CIRI2", 
            "gCSI_CIRCexplorer2", "CCLE_CIRCexplorer2", "GDSC_CIRCexplorer2",
            "gCSI_circRNA_finder", "CCLE_circRNA_finder", "GDSC_circRNA_finder", 
            "gCSI_find_circ", "CCLE_find_circ", "GDSC_find_circ")
set_order2 = c("gCSI_CIRI2", "gCSI_CIRCexplorer2", "gCSI_circRNA_finder", "gCSI_find_circ",
            "CCLE_CIRCexplorer2", "CCLE_CIRI2", "CCLE_circRNA_finder", "CCLE_find_circ",
            "GDSC_CIRI2", "GDSC_CIRCexplorer2", "GDSC_circRNA_finder", "GDSC_find_circ")

# plot upset plots
pdf("../results/figures/figure4/upset_by_pipeline_ge.pdf", width=10, height=5)
plot_upset(toPlot, set_order = set_order1)
dev.off()

pdf("../results/figures/figure4/upset_by_pset_ge.pdf", width=10, height=5)
plot_upset(toPlot, set_order = set_order2)
dev.off()


############################################################
# Figure 4: Plot transcript exp counts per pipeline
############################################################

# create data frame of counts for plotting
df <- data.frame(Count = c(sum(ciri_gcsi_ge), sum(ciri_ccle_ge), sum(ciri_gdsc_ge),
                           sum(circ_gcsi_ge), sum(circ_ccle_ge), sum(circ_gdsc_ge), 
                           sum(cfnd_gcsi_ge), sum(cfnd_ccle_ge), sum(cfnd_gdsc_ge), 
                           sum(fcrc_gcsi_ge), sum(fcrc_ccle_ge), sum(fcrc_gdsc_ge)),
                PSet = c(rep(c("gCSI", "CCLE", "GDSC2"), 4)),
                Pipeline = c(rep("CIRI2", 3), rep("CIRCexplorer2", 3), rep("circRNA_finder", 3), rep("find_circ", 3)))

df$Pipeline <- factor(df$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
df$PSet <- factor(df$PSet, levels = c("gCSI", "CCLE", "GDSC2"))

# plot bar plot of counts
png("../results/figures/figure4/counts_ge.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(df, aes(x = Pipeline, y = Count, fill = PSet)) + geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_fill_manual(values=pal1, limits=c("gCSI", "CCLE", "GDSC2")) + 
  scale_y_continuous(limits = c(0, 240000), expand=c(0,0))  + theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.4, 'cm'))
dev.off()


############################################################
# Figure 4: Plot distribution of transcript exp
############################################################

# create data frame of expression values
print("making expression")
Expression = c(as.vector(unlist(ciri_gcsi_ge)), as.vector(unlist(ciri_ccle_ge)), as.vector(unlist(ciri_gdsc_ge)),
                as.vector(unlist(circ_gcsi_ge)), as.vector(unlist(circ_ccle_ge)), as.vector(unlist(circ_gdsc_ge)),
                as.vector(unlist(cfnd_gcsi_ge)), as.vector(unlist(cfnd_ccle_ge)), as.vector(unlist(cfnd_gdsc_ge)),
                as.vector(unlist(fcrc_gcsi_ge)), as.vector(unlist(fcrc_ccle_ge)), as.vector(unlist(fcrc_gdsc_ge)))
print("making label")
Label = c(rep("gCSI-CIRI2", prod(dim(ciri_gcsi_ge))), rep("CCLE-CIRI2", prod(dim(ciri_ccle_ge))), rep("GDSC2-CIRI2", prod(dim(ciri_gdsc_ge))),
        rep("gCSI-CIRCexplorer2", prod(dim(circ_gcsi_ge))), rep("CCLE-CIRCexplorer2", prod(dim(circ_ccle_ge))), rep("GDSC2-CIRCexplorer2", prod(dim(circ_gdsc_ge))),
        rep("gCSI-circRNA_finder", prod(dim(cfnd_gcsi_ge))), rep("CCLE-circRNA_finder", prod(dim(cfnd_ccle_ge))), rep("GDSC2-circRNA_finder", prod(dim(cfnd_gdsc_ge))),
        rep("gCSI-find_circ", prod(dim(fcrc_gcsi_ge))), rep("CCLE-find_circ", prod(dim(fcrc_ccle_ge))), rep("GDSC2-find_circ", prod(dim(fcrc_gdsc_ge))))
                
toPlot <- data.frame(Expression, Label)
dim(toPlot)
toPlot$Pipeline <- gsub(".*-", "", toPlot$Label)
toPlot$PSet <- gsub("-.*", "", toPlot$Label)

toPlot$Pipeline <- factor(toPlot$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
toPlot$PSet <- factor(toPlot$PSet, levels = c("gCSI", "CCLE", "GDSC2"))


max = max(toPlot$Expression)

png("../results/figures/figure4/count_dist_all_ge.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure4/count_dist_pipeline_ge.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure4/count_dist_pset_ge.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

