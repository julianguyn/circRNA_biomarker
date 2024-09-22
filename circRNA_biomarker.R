# load in packages
suppressMessages(library(data.table))
#suppressMessages(library(PharmacoGx))
suppressMessages(library(survcomp))

# read in circRNA matrices
suppressWarnings(ciri_gcsi <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gcsi_counts.tsv", data.table = F))
suppressWarnings(ciri_gdsc <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gdsc_counts.tsv", data.table = F))
suppressWarnings(ciri_ccle <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_ccle_counts.tsv", data.table = F))

suppressWarnings(circ_gcsi <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F))
suppressWarnings(circ_gdsc <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F))
suppressWarnings(circ_ccle <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_ccle_counts.tsv", data.table = F))

suppressWarnings(cfnd_gcsi <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F))
suppressWarnings(cfnd_gdsc <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F))
suppressWarnings(cfnd_ccle <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_ccle_counts.tsv", data.table = F))

suppressWarnings(fcrc_gcsi <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gcsi_counts.tsv", data.table = F))
suppressWarnings(fcrc_gdsc <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gdsc_counts.tsv", data.table = F))
suppressWarnings(fcrc_ccle <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_ccle_counts.tsv", data.table = F))


# TODO: load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

# load in metadata
load("../results/data/tissue-metadata.RData")

# match drug sensitivity cell line names to cellid names
# TODO: manually check for cell line overlap (skipped) -- gcsi 334 in, 75 not - ctrp: 821 in, 66 not, gdsc: 329 in, 480 not
gcsi_sen <- gcsi_sen[,colnames(gcsi_sen) %in% gcsi$cellid]
ctrp_sen <- ctrp_sen[,colnames(ctrp_sen) %in% ccle$cellid]
gdsc_sen <- gdsc_sen[,colnames(gdsc_sen) %in% gdsc$cellid]

# TODO: run on all drugs from cRic
drugs <- intersect(intersect(rownames(gcsi_sen), rownames(ctrp_sen)), rownames(gdsc_sen))
# "Bortezomib", "Crizotinib", "Docetaxel", "Erlotinib", "Pictilisib", "Gemcitabine", "Lapatinib", "Entinostat", "Paclitaxel", "Sirolimus", "Vorinostat" 
gcsi_sen <- gcsi_sen[rownames(gcsi_sen) %in% drugs,]
ctrp_sen <- ctrp_sen[rownames(ctrp_sen) %in% drugs,]
gdsc_sen <- gdsc_sen[rownames(gdsc_sen) %in% drugs,]


## Compute Concordance Index
computeCI <- function(circ_counts, sensitivity_data, pset_label, pipeline_label){

    print(paste0("-------- STARTING: ", pset_label, "-", pipeline_label))

    # filter to only keep common cell lines
    common_samples <- intersect(circ_counts$sample, colnames(sensitivity_data))
    circ_counts <- circ_counts[circ_counts$sample %in% common_samples,]
    sensitivity_data <- sensitivity_data[,colnames(sensitivity_data) %in% common_samples]

    # order dataframes
    circ_counts <- circ_counts[order(circ_counts$sample),]
    sensitivity_data <- sensitivity_data[,order(colnames(sensitivity_data))]

    # format circRNA counts matrix
    samples <- circ_counts$sample
    circ_counts$sample <- NULL
    circ_counts <- as.data.frame(t(circ_counts))
    colnames(circ_counts) <- samples

    # remove circRNA transcripts with 0 expression across common samples
    circ_counts[is.na(circ_counts)] <- 0                                                            ## TODO: remove - check where NAs are from
    circ_counts <- circ_counts[rowSums(circ_counts[])>0,]

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(circ_counts) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("circRNA", "drug", "ci", "pvalue", "se", "upper", "lower")
    combinations$circRNA <- rep(rownames(circ_counts), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(circ_counts))

    print(nrow(combinations))
  
    # compute concordance index
    for (i in 1:nrow(combinations)){
        if (i %% 10000 == 0) {print(paste0(i, " out of ", nrow(combinations), " complete"))}

        ci <- survcomp::concordance.index(as.numeric(sensitivity_data[combinations[,2][i],]), 
                                            # row from sensitivity_data with the drug for all samples
                                            surv.time = as.numeric(unlist(-circ_counts[combinations[,1][i],])), 
                                            # row of circRNA transcript i in circRNA expression matrix with cell lines as columns
                                            surv.event = rep(1,length(sensitivity_data)), 
                                            # df of all drugs as rows with all samples as columns
                                            outx = TRUE, method="noether", na.rm = TRUE)

        combinations$pvalue[i] <- ci$p.value
        combinations$ci[i] <- ci$c.index
        combinations$se[i] <- ci$se
        combinations$upper[i] <- ci$upper
        combinations$lower[i] <- ci$lower

        
    }

    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    combinations$FDR <- p.adjust(combinations$pvalue, method = "BH", n = length(combinations$pvalue))
    combinations$FDRsig <- ifelse(combinations$FDR < 0.05, TRUE, FALSE)

    
    # format dataframe for plotting (order by ci and add variable rank)
    combinations <- combinations[order(combinations$ci),]
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$circRNA, "-", combinations$drug)
    combinations$pset <- pset_label
    combinations$pipeline <- pipeline_label
    
    return(combinations)
}

gcsi_ci_ciri <- computeCI(ciri_gcsi, gcsi_sen, "gCSI", "CIRI2")
ccle_ci_ciri <- computeCI(ciri_ccle, ccle_sen, "CCLE", "CIRI2")
gdsc_ci_ciri <- computeCI(ciri_gdsc, gdsc_sen, "GDSC2", "CIRI2")

save(gcsi_ci_ciri, ccle_ci_ciri, gdsc_ci_ciri, file = "../reuslts/data/ci/ciri_ci.RData")

gcsi_ci_circ <- computeCI(circ_gcsi, gcsi_sen, "gCSI", "CIRCexplorer2")
ccle_ci_circ <- computeCI(circ_ccle, ccle_sen, "CCLE", "CIRCexplorer2")
gdsc_ci_circ <- computeCI(circ_gdsc, gdsc_sen, "GDSC2", "CIRCexplorer2")

save(gcsi_ci_circ, ccle_ci_circ, gdsc_ci_circ, file = "../reuslts/data/ci/circ_ci.RData")

gcsi_ci_cfnd <- computeCI(cfnd_gcsi, gcsi_sen, "gCSI", "circRNA_finder")
ccle_ci_cfnd <- computeCI(cfnd_ccle, ccle_sen, "CCLE", "circRNA_finder")
gdsc_ci_cfnd <- computeCI(cfnd_gdsc, gdsc_sen, "GDSC2", "circRNA_finder")

save(gcsi_ci_cfnd, ccle_ci_cfnd, gdsc_ci_cfnd, file = "../reuslts/data/ci/cfnd_ci.RData")

gcsi_ci_fcrc <- computeCI(fcrc_gcsi, gcsi_sen, "gCSI", "find_circ")
ccle_ci_fcrc <- computeCI(fcrc_ccle, ccle_sen, "CCLE", "find_circ")
gdsc_ci_fcrc <- computeCI(fcrc_gdsc, gdsc_sen, "GDSC2", "find_circ")

save(gcsi_ci_fcrc, ccle_ci_fcrc, gdsc_ci_fcrc, file = "../reuslts/data/ci/fcrc_ci.RData")



