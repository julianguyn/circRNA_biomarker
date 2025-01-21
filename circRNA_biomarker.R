# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(survcomp)
    library(ggplot2)
    library(ggpubr)
    library(ggvenn)
    library(meta)
    library(ComplexHeatmap)
})

############################################################
# Load in data 
############################################################

# load circRNA expression data
path <- "../data/processed_cellline/GE_common_samples/"

ciri_gcsi <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F)                          
ciri_gdsc <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)                          
ciri_ccle <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)                          

circ_gcsi <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)                  
circ_gdsc <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)                  
circ_ccle <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F)                  

cfnd_gcsi <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)                 
cfnd_gdsc <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)                 
cfnd_ccle <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F)                 

fcrc_gcsi <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)                      
fcrc_gdsc <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)                      
fcrc_ccle <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F)                        

# load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

# This is for all_samples if we go down that route
# load in metadata
#load("../results/data/tissue-metadata.RData")


############################################################
# Subset data to samples and drugs of interest
############################################################

# drugs of interest
drugs <- intersect(intersect(rownames(gcsi_sen), rownames(ctrp_sen)), rownames(gdsc_sen))
# "Bortezomib", "Crizotinib", "Docetaxel", "Erlotinib", "Pictilisib", "Gemcitabine", "Lapatinib", "Entinostat", "Paclitaxel", "Sirolimus", "Vorinostat" 

# common samples
samples <- ciri_gcsi$sample

# keep common samples and drugs of interest
gcsi_sen <- gcsi_sen[rownames(gcsi_sen) %in% drugs, colnames(gcsi_sen) %in% samples]
ctrp_sen <- ctrp_sen[rownames(ctrp_sen) %in% drugs, colnames(ctrp_sen) %in% samples]
gdsc_sen <- gdsc_sen[rownames(gdsc_sen) %in% drugs, colnames(gdsc_sen) %in% samples]

# This is for all_samples if we go down that route
# match drug sensitivity cell line names to cellid names
# TODO: manually check for cell line overlap (skipped) -- gcsi 334 in, 75 not - ctrp: 821 in, 66 not, gdsc: 329 in, 480 not
#gcsi_sen <- gcsi_sen[,colnames(gcsi_sen) %in% gcsi$cellid]
#ctrp_sen <- ctrp_sen[,colnames(ctrp_sen) %in% ccle$cellid]
#gdsc_sen <- gdsc_sen[,colnames(gdsc_sen) %in% gdsc$cellid]



############################################################
# Compute concordance index 
############################################################

# function to compute concordance index
computeCI <- function(circ_counts, sensitivity_data, pset_label, pipeline_label){ 

        print(paste0("-------- STARTING: ", pset_label, "-", pipeline_label))

        # format circ_counts dataframe
        rownames(circ_counts) <- circ_counts$sample
        circ_counts$sample <- NULL
        circ_counts <- t(circ_counts) |> as.data.frame()
        
        # keep and order common cell lines
        common_cells <- intersect(colnames(circ_counts), colnames(sensitivity_data))
        circ_counts <- circ_counts[,common_cells]
        sensitivity_data <- sensitivity_data[,common_cells]

        # keep only circRNA transcripts that have expression in at least 3 cell lines
        circ_counts[is.na(circ_counts)] <- 0                                                            ## TODO: remove - check where NAs are from
        circ_counts <- circ_counts[rowSums(circ_counts[])>2,]

        # create data frame to hold results
        combinations <- expand.grid(circRNA = rownames(circ_counts), drug = rownames(sensitivity_data))
        num_combinations <- nrow(combinations)

        if (nrow(circ_counts) > 0) {
                # initialize columns to store univariable results
                combinations$ci <- NA
                combinations$pvalue <- NA
                combinations$se <- NA
                combinations$upper <- NA
                combinations$lower <- NA

                # extract drug and feature data into lists
                sen_list <- lapply(rownames(sensitivity_data), function(drug) as.numeric(sensitivity_data[drug, ])) # each is drug response for one drug
                feat_list <- lapply(rownames(circ_counts), function(circRNA) as.numeric(circ_counts[circRNA, ])) # each is expression for each transcript

                # compute concordance index for each feature-drug pair
                for (i in 1:num_combinations) {
                        if (i %% 10000 == 0) { print(paste0(i, " out of ", num_combinations, " complete"))}
                        
                        # compute concordance index
                        ci <- survcomp::concordance.index(
                        sen_list[[which(rownames(sensitivity_data) == combinations$drug[i])]],                  # sensitivity data for all samples
                        surv.time = feat_list[[which(rownames(circ_counts) == combinations$circRNA[i])]],          # feature values for all samples
                        surv.event = rep(1, ncol(circ_counts)),
                        outx = TRUE,
                        method = "noether",
                        na.rm = TRUE
                        )

                        combinations$ci[i] <- ci$c.index
                        combinations$pvalue[i] <- ci$p.value
                        combinations$se[i] <- ci$se
                        combinations$upper[i] <- ci$upper
                        combinations$lower[i] <- ci$lower
                }

                # filtering and multiple test correction
                #combinations <- combinations[complete.cases(combinations$pvalue),]
                #combinations$FDR <- p.adjust(combinations$pvalue, method = "BH")
                #combinations$FDRsig <- combinations$FDR < 0.05

                # format dataframe for plotting (order by ci and add variable rank)
                combinations <- combinations[order(combinations$ci),]
                combinations$rank <- seq_len(nrow(combinations))
                combinations$pairs <- paste0(combinations$circRNA, "-", combinations$drug)
                combinations$pset <- pset_label
                combinations$pipeline <- pipeline_label
        }

        print(dim(combinations))

        return(combinations)
}

gcsi_ci_ciri <- computeCI(ciri_gcsi, gcsi_sen, "gCSI", "CIRI2")
ccle_ci_ciri <- computeCI(ciri_ccle, ccle_sen, "CCLE", "CIRI2")
gdsc_ci_ciri <- computeCI(ciri_gdsc, gdsc_sen, "GDSC2", "CIRI2")

gcsi_ci_circ <- computeCI(circ_gcsi, gcsi_sen, "gCSI", "CIRCexplorer2")
ccle_ci_circ <- computeCI(circ_ccle, ccle_sen, "CCLE", "CIRCexplorer2")
gdsc_ci_circ <- computeCI(circ_gdsc, gdsc_sen, "GDSC2", "CIRCexplorer2")

gcsi_ci_cfnd <- computeCI(cfnd_gcsi, gcsi_sen, "gCSI", "circRNA_finder")
ccle_ci_cfnd <- computeCI(cfnd_ccle, ccle_sen, "CCLE", "circRNA_finder")
gdsc_ci_cfnd <- computeCI(cfnd_gdsc, gdsc_sen, "GDSC2", "circRNA_finder")

gcsi_ci_fcrc <- computeCI(fcrc_gcsi, gcsi_sen, "gCSI", "find_circ")
ccle_ci_fcrc <- computeCI(fcrc_ccle, ccle_sen, "CCLE", "find_circ")
gdsc_ci_fcrc <- computeCI(fcrc_gdsc, gdsc_sen, "GDSC2", "find_circ")


############################################################
# Save results
############################################################

save(gcsi_ci_ciri, ccle_ci_ciri, gdsc_ci_ciri, file = "../results/data/ci/GE_ciri_ci.RData")
save(gcsi_ci_circ, ccle_ci_circ, gdsc_ci_circ, file = "../results/data/ci/GE_circ_ci.RData")
save(gcsi_ci_cfnd, ccle_ci_cfnd, gdsc_ci_cfnd, file = "../results/data/ci/GE_cfnd_ci.RData")
save(gcsi_ci_fcrc, ccle_ci_fcrc, gdsc_ci_fcrc, file = "../results/data/ci/GE_fcrc_ci.RData")



############################################################
# Load data (separate script on h4h)
############################################################

load("../reuslts/data/ci/GE_ciri_ci1_1.RData")
load("../reuslts/data/ci/GE_ciri_ci1_2.RData")
load("../reuslts/data/ci/GE_ciri_ci1_3.RData")

load("../reuslts/data/ci/GE_circ_ci2_1.RData")
load("../reuslts/data/ci/GE_circ_ci2_2.RData")
load("../reuslts/data/ci/GE_circ_ci2_3.RData")

load("../reuslts/data/ci/GE_cfnd_ci3_1.RData")
load("../reuslts/data/ci/GE_cfnd_ci3_2.RData")
load("../reuslts/data/ci/GE_cfnd_ci3_3.RData")

load("../reuslts/data/ci/GE_fcrc_ci4_1.RData")
load("../reuslts/data/ci/GE_fcrc_ci4_2.RData")
load("../reuslts/data/ci/GE_fcrc_ci4_3.RData")



##### Compare number of significant predictions per dataset per pipeline
ciri_transcripts <- list(gCSI = gcsi_ci_ciri[which(gcsi_ci_ciri$FDR < 0.05 & (gcsi_ci_ciri$ci > 0.65 | gcsi_ci_ciri$ci < 0.35)),]$pairs, 
                         CTRP = ccle_ci_ciri[which(ccle_ci_ciri$FDR < 0.05 & (ccle_ci_ciri$ci > 0.65 | ccle_ci_ciri$ci < 0.35)),]$pairs, 
                         GDSC2 = gdsc_ci_ciri[which(gdsc_ci_ciri$FDR < 0.05 & (gdsc_ci_ciri$ci > 0.65 | gdsc_ci_ciri$ci < 0.35)),]$pairs)
circ_transcripts <- list(gCSI = gcsi_ci_circ[which(gcsi_ci_circ$FDR < 0.05 & (gcsi_ci_circ$ci > 0.65 | gcsi_ci_circ$ci < 0.35)),]$pairs, 
                         CTRP = ccle_ci_circ[which(ccle_ci_circ$FDR < 0.05 & (ccle_ci_circ$ci > 0.65 | ccle_ci_circ$ci < 0.35)),]$pairs, 
                         GDSC2 = gdsc_ci_circ[which(gdsc_ci_circ$FDR < 0.05 & (gdsc_ci_circ$ci > 0.65 | gdsc_ci_circ$ci < 0.35)),]$pairs)
cfnd_transcripts <- list(gCSI = gcsi_ci_cfnd[which(gcsi_ci_cfnd$FDR < 0.05 & (gcsi_ci_cfnd$ci > 0.65 | gcsi_ci_cfnd$ci < 0.35)),]$pairs, 
                         CTRP = ccle_ci_cfnd[which(ccle_ci_cfnd$FDR < 0.05 & (ccle_ci_cfnd$ci > 0.65 | ccle_ci_cfnd$ci < 0.35)),]$pairs, 
                         GDSC2 = gdsc_ci_cfnd[which(gdsc_ci_cfnd$FDR < 0.05 & (gdsc_ci_cfnd$ci > 0.65 | gdsc_ci_cfnd$ci < 0.35)),]$pairs)
fcrc_transcripts <- list(gCSI = gcsi_ci_fcrc[which(gcsi_ci_fcrc$FDR < 0.05 & (gcsi_ci_fcrc$ci > 0.65 | gcsi_ci_fcrc$ci < 0.35)),]$pairs, 
                         CTRP = ccle_ci_fcrc[which(ccle_ci_fcrc$FDR < 0.05 & (ccle_ci_fcrc$ci > 0.65 | ccle_ci_fcrc$ci < 0.35)),]$pairs, 
                         GDSC2 = gdsc_ci_fcrc[which(gdsc_ci_fcrc$FDR < 0.05 & (gdsc_ci_fcrc$ci > 0.65 | gdsc_ci_fcrc$ci < 0.35)),]$pairs)

# plot venn diagram
p1 <- ggvenn(ciri_transcripts, 
        fill_color = c("#51C7AD", "#392C57", "#3670A0"), stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "CIRI2")
p2 <- ggvenn(circ_transcripts, 
        fill_color = c("#51C7AD", "#392C57", "#3670A0"), stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "\nCIRCexplorer2")
p3 <- ggvenn(cfnd_transcripts, 
        fill_color = c("#51C7AD", "#392C57", "#3670A0"), stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "\ncircRNA_finder")
p4 <- ggvenn(fcrc_transcripts, 
        fill_color = c("#51C7AD", "#392C57", "#3670A0"), stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "\nfind_circ")

png("../results/figures/figure3/venndiagram_per_pipeline.png", width=150, height=500, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4, common.legend = FALSE)
dev.off()

##################
### Upset plot ###
##################

set.seed(123)
toPlot <- make_comb_mat(list(
            gCSI_CIRI2 = gcsi_ci_ciri[which(gcsi_ci_ciri$FDR < 0.05 & (gcsi_ci_ciri$ci > 0.65 | gcsi_ci_ciri$ci < 0.35)),]$pairs,
            CTRP_CIRI2 = ccle_ci_ciri[which(ccle_ci_ciri$FDR < 0.05 & (ccle_ci_ciri$ci > 0.65 | ccle_ci_ciri$ci < 0.35)),]$pairs,
            GDSC_CIRI2 = gdsc_ci_ciri[which(gdsc_ci_ciri$FDR < 0.05 & (gdsc_ci_ciri$ci > 0.65 | gdsc_ci_ciri$ci < 0.35)),]$pairs,
            gCSI_CIRCexplorer2 = gcsi_ci_circ[which(gcsi_ci_circ$FDR < 0.05 & (gcsi_ci_circ$ci > 0.65 | gcsi_ci_circ$ci < 0.35)),]$pairs,
            CTRP_CIRCexplorer2 = ccle_ci_circ[which(ccle_ci_circ$FDR < 0.05 & (ccle_ci_circ$ci > 0.65 | ccle_ci_circ$ci < 0.35)),]$pairs,
            GDSC_CIRCexplorer2 = gdsc_ci_circ[which(gdsc_ci_circ$FDR < 0.05 & (gdsc_ci_circ$ci > 0.65 | gdsc_ci_circ$ci < 0.35)),]$pairs,
            gCSI_circRNA_finder = gcsi_ci_cfnd[which(gcsi_ci_cfnd$FDR < 0.05 & (gcsi_ci_cfnd$ci > 0.65 | gcsi_ci_cfnd$ci < 0.35)),]$pairs,
            CTRP_circRNA_finder = ccle_ci_cfnd[which(ccle_ci_cfnd$FDR < 0.05 & (ccle_ci_cfnd$ci > 0.65 | ccle_ci_cfnd$ci < 0.35)),]$pairs,
            GDSC_circRNA_finder = gdsc_ci_cfnd[which(gdsc_ci_cfnd$FDR < 0.05 & (gdsc_ci_cfnd$ci > 0.65 | gdsc_ci_cfnd$ci < 0.35)),]$pairs,
            gCSI_find_circ = gcsi_ci_fcrc[which(gcsi_ci_fcrc$FDR < 0.05 & (gcsi_ci_fcrc$ci > 0.65 | gcsi_ci_fcrc$ci < 0.35)),]$pairs,
            CTRP_find_circ = ccle_ci_fcrc[which(ccle_ci_fcrc$FDR < 0.05 & (ccle_ci_fcrc$ci > 0.65 | ccle_ci_fcrc$ci < 0.35)),]$pairs,
            GDSC_find_circ = gdsc_ci_fcrc[which(gdsc_ci_fcrc$FDR < 0.05 & (gdsc_ci_fcrc$ci > 0.65 | gdsc_ci_fcrc$ci < 0.35)),]$pairs))

# remove combinations of less than 2 pairs
toPlot <- toPlot[comb_size(toPlot) >= 2]

# upset plot
pdf("../results/figures/figure3/upset_by_pipeline.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "CTRP_CIRI2", "GDSC_CIRI2", "gCSI_CIRCexplorer2", "CTRP_CIRCexplorer2", "GDSC_CIRCexplorer2",
                            "gCSI_circRNA_finder", "CTRP_circRNA_finder", "GDSC_circRNA_finder", "gCSI_find_circ", "CTRP_find_circ", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()

pdf("../results/figures/figure3/upset_by_pset.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "gCSI_CIRCexplorer2", "gCSI_circRNA_finder", "gCSI_find_circ",
                            "CTRP_CIRI2", "CTRP_CIRCexplorer2", "CTRP_circRNA_finder", "CTRP_find_circ",
                            "GDSC_CIRI2", "GDSC_CIRCexplorer2", "GDSC_circRNA_finder", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()

