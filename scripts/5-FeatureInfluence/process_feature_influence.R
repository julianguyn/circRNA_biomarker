# load libraries
suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(tidyverse)
    library(ggh4x)
    library(RColorBrewer)
    library(PharmacoGx)
})

set.seed(200)


############################################################
# Load in data
############################################################

# load in stability data
load("../results/data/temp/gene_isoform_stability.RData")
load("../results/data/temp/circ_stability2.RData")

# load in circRNA features
load("../results/data/circ_stability_gcsi_features2.RData")
load("../results/data/circ_stability_ccle_features2.RData")
load("../results/data/circ_stability_gdsc_features2.RData")


############################################################
# Format circRNA dataframes
############################################################

# function to format circRNA stability dataframes
c_stability <- function(stability, gcsi, ccle, gdsc) {

  # reformat stability dataframe
  stability <- reshape2::dcast(stability, transcript ~ pair, value.var = "si")
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
  feats <- data.frame(transcript = c(rownames(gcsi), rownames(ccle), rownames(gdsc)),
                      gc = c(gcsi$GC, ccle$GC, gdsc$GC),
                      nexons = c(gcsi$NExons, ccle$NExons, gdsc$NExons),
                      length = c(gcsi$Length, ccle$Length, gdsc$Length))
  feats <- unique(feats)
  feats <- feats[match(stability$transcript, feats$transcript),]
  
  stability$gc <- feats$gc
  stability$n_exon <- feats$nexons
  stability$length <- feats$length
  stability$chr <- gsub("\\..*", "", stability$transcript)

  return(stability)
} 

ciri_stability <- c_stability(ciri_stability, ciri_gcsi_ft, ciri_ccle_ft, ciri_gdsc_ft)
circ_stability <- c_stability(circ_stability, circ_gcsi_ft, circ_ccle_ft, circ_gdsc_ft)
cfnd_stability <- c_stability(cfnd_stability, cfnd_gcsi_ft, cfnd_ccle_ft, cfnd_gdsc_ft)
fcrc_stability <- c_stability(fcrc_stability, fcrc_gcsi_ft, fcrc_ccle_ft, fcrc_gdsc_ft)

# save files
write.csv(ciri_stability, file = "../results/data/temp/ciri_stability.csv", quote = F, row.names = F)
write.csv(circ_stability, file = "../results/data/temp/circ_stability.csv", quote = F, row.names = F)
write.csv(cfnd_stability, file = "../results/data/temp/cfnd_stability.csv", quote = F, row.names = F)
write.csv(fcrc_stability, file = "../results/data/temp/fcrc_stability.csv", quote = F, row.names = F)