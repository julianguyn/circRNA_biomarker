# load libraries
suppressPackageStartupMessages({
    library(caret)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
    library(dplyr)
})

set.seed(200)  


############################################################
# Load in data
############################################################

# load in stability data
load("../results/data/temp/gene_isoform_stability.RData")
load("../results/data/temp/circ_stability.RData")

# load in circRNA features
load("../results/data/circ_stability_gcsi_features.RData")
load("../results/data/circ_stability_ccle_features.RData")
load("../results/data/circ_stability_gdsc_features.RData")


############################################################
# Format circRNA dataframes
############################################################

# function to format circRNA stability dataframes
c_stability <- function(stability, gcsi, ccle, gdsc) {
  colnames(stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman", "circ_id")
  stability$gcsi_median <- gcsi$MedianExp
  stability$ccle_median <- ccle$MedianExp
  stability$gdsc_median <- gdsc$MedianExp
  stability$gc <- gcsi$GC
  stability$n_exon <- gcsi$NExons
  stability$length <- gcsi$Length
  return(stability)
} 

ciri_stability <- c_stability(ciri_stability, ciri_gcsi_ft, ciri_ccle_ft, ciri_gdsc_ft)
circ_stability <- c_stability(circ_stability, circ_gcsi_ft, circ_ccle_ft, circ_gdsc_ft)
cfnd_stability <- c_stability(cfnd_stability, cfnd_gcsi_ft, cfnd_ccle_ft, cfnd_gdsc_ft)
fcrc_stability <- c_stability(fcrc_stability, fcrc_gcsi_ft, fcrc_ccle_ft, fcrc_gdsc_ft)


############################################################
# Specify quantile based on Kuncheva Index
############################################################

gene_k <- 0.25
isof_k <- 0.15
ciri_k <- 0.2
circ_k <- 0.1
cfnd_k <- 0.15
fcrc_k <- 0.1

############################################################
# Create dataset pairs 
############################################################

# function to split stability dataframe for each dataset pair
stability_pair <- function(stability, pair, kuncheva) {

  n <- nrow(stability)
  
  if (pair == "gcsi_ccle") {
    stability <- stability[,c('gdsc_median','gc','n_exon','length','gcsi_ccle_spearman')]
    stable <- stability[order(-stability$gcsi_ccle_spearman),]
    unstable <- stability[order(stability$gcsi_ccle_spearman),]
  }
  if (pair == "gcsi_gdsc") {
    stability <- stability[,c('ccle_median','gc','n_exon','length','gcsi_gdsc_spearman')]
    stable <- stability[order(-stability$gcsi_gdsc_spearman),]
    unstable <- stability[order(stability$gcsi_gdsc_spearman),]
  }
  if (pair == "gdsc_ccle") {
    stability <- stability[,c('gcsi_median','gc','n_exon','length','gdsc_ccle_spearman')]
    stable <- stability[order(-stability$gdsc_ccle_spearman),]
    unstable <- stability[order(stability$gdsc_ccle_spearman),]
  }

  # subset based on kuncheva index
  stable <- stable[1:(kuncheva*n),]
  unstable <- unstable[1:(kuncheva*n),]
  stability <- rbind(stable,unstable)

  # shuffle data
  stability <- stability[sample(1:nrow(stability)),]

  return(stability)
}

gene_stability_gcsi_ccle <- stability_pair(gene_stability, "gcsi_ccle", gene_k)
gene_stability_gcsi_gdsc <- stability_pair(gene_stability, "gcsi_gdsc", gene_k)
gene_stability_gdsc_ccle <- stability_pair(gene_stability, "gdsc_ccle", gene_k)

transcript_stability_gcsi_ccle <- stability_pair(transcript_stability, "gcsi_ccle", isof_k)
transcript_stability_gcsi_gdsc <- stability_pair(transcript_stability, "gcsi_gdsc", isof_k)
transcript_stability_gdsc_ccle <- stability_pair(transcript_stability, "gdsc_ccle", isof_k)

ciri_stability_gcsi_ccle <- stability_pair(ciri_stability, "gcsi_ccle", ciri_k) #12
ciri_stability_gcsi_gdsc <- stability_pair(ciri_stability, "gcsi_gdsc", ciri_k)
ciri_stability_gdsc_ccle <- stability_pair(ciri_stability, "gdsc_ccle", ciri_k)

circ_stability_gcsi_ccle <- stability_pair(circ_stability, "gcsi_ccle", circ_k) #12
circ_stability_gcsi_gdsc <- stability_pair(circ_stability, "gcsi_gdsc", circ_k)
circ_stability_gdsc_ccle <- stability_pair(circ_stability, "gdsc_ccle", circ_k)

cfnd_stability_gcsi_ccle <- stability_pair(cfnd_stability, "gcsi_ccle", cfnd_k) #24
cfnd_stability_gcsi_gdsc <- stability_pair(cfnd_stability, "gcsi_gdsc", cfnd_k)
cfnd_stability_gdsc_ccle <- stability_pair(cfnd_stability, "gdsc_ccle", cfnd_k)

fcrc_stability_gcsi_ccle <- stability_pair(fcrc_stability, "gcsi_ccle", fcrc_k) #100
fcrc_stability_gcsi_gdsc <- stability_pair(fcrc_stability, "gcsi_gdsc", fcrc_k)
fcrc_stability_gdsc_ccle <- stability_pair(fcrc_stability, "gdsc_ccle", fcrc_k)



############################################################
# Permutation feature importance
############################################################

# function to run Linear Regression
runLinearReg <- function(x, label) {
  
  # Create dataframe to store feature scores
  feature_score <- data.frame(matrix(ncol=4, nrow = 0))
  colnames(feature_score) <- c("feature", "baseline", "permuted", "dataset")
  
  # Change colnames of inputted dataframe
  colnames(x) <- c("median", colnames(x)[2:4], "spearman")
  x <- na.omit(x)
  
  #10-fold CV
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10) 
  
  # Loop through each feature
  for (i in 1:4) {
    
    feature <- colnames(x)[i]
    
    # get baseline score
    gbmFit <- train(spearman ~ ., data = x, 
                    method = "lm", 
                    trControl = fitControl,
                    verbose = FALSE)
    #Save MSE score
    baseline_score <- (gbmFit$results$RMSE)^2
    
    # Permute each feature 20,000 times
    for (k in 1:20000) { 
      
      # shuffle feature column
      shuffle_feat <- sample(x[,i])
      x[,i] <- shuffle_feat
      
      gbmFit <- train(spearman ~ ., data = x, 
                      method = "lm", 
                      trControl = fitControl,
                      verbose = FALSE)
      
      # Record permuted MSE score
      permutation_score <- (gbmFit$results$RMSE)^2
      
      score <- data.frame(feature = feature, baseline = baseline_score, permuted = permutation_score, dataset = label)
      feature_score <- rbind(feature_score, score)
    }
  }
  return(feature_score)
}


gene_stability_gcsi_ccle <- runLinearReg(gene_stability_gcsi_ccle, "gcsi_ccle")
gene_stability_gcsi_gdsc <- runLinearReg(gene_stability_gcsi_gdsc, "gcsi_gdsc")
gene_stability_gdsc_ccle <- runLinearReg(gene_stability_gdsc_ccle, "gdsc_ccle")

transcript_stability_gcsi_ccle <- runLinearReg(transcript_stability_gcsi_ccle, "gcsi_ccle")
transcript_stability_gcsi_gdsc <- runLinearReg(transcript_stability_gcsi_gdsc, "gcsi_gdsc")
transcript_stability_gdsc_ccle <- runLinearReg(transcript_stability_gdsc_ccle, "gdsc_ccle")

ciri_stability_gcsi_ccle <- runLinearReg(ciri_stability_gcsi_ccle, "gcsi_ccle")
ciri_stability_gcsi_gdsc <- runLinearReg(ciri_stability_gcsi_gdsc, "gcsi_gdsc")
ciri_stability_gdsc_ccle <- runLinearReg(ciri_stability_gdsc_ccle, "gdsc_ccle")

circ_stability_gcsi_ccle <- runLinearReg(circ_stability_gcsi_ccle, "gcsi_ccle")
circ_stability_gcsi_gdsc <- runLinearReg(circ_stability_gcsi_gdsc, "gcsi_gdsc")
circ_stability_gdsc_ccle <- runLinearReg(circ_stability_gdsc_ccle, "gdsc_ccle")

cfnd_stability_gcsi_ccle <- runLinearReg(cfnd_stability_gcsi_ccle, "gcsi_ccle")
cfnd_stability_gcsi_gdsc <- runLinearReg(cfnd_stability_gcsi_gdsc, "gcsi_gdsc")
cfnd_stability_gdsc_ccle <- runLinearReg(cfnd_stability_gdsc_ccle, "gdsc_ccle")

fcrc_stability_gcsi_ccle <- runLinearReg(fcrc_stability_gcsi_ccle, "gcsi_ccle")
fcrc_stability_gcsi_gdsc <- runLinearReg(fcrc_stability_gcsi_gdsc, "gcsi_gdsc")
fcrc_stability_gdsc_ccle <- runLinearReg(fcrc_stability_gdsc_ccle, "gdsc_ccle")

############################################################
# Combine results for each transcript type
############################################################

gene_features <- rbind(gene_stability_gcsi_ccle, gene_stability_gcsi_gdsc, gene_stability_ggdsc_ccle)
transcript_features <- rbind(transcript_stability_gcsi_ccle, transcript_stability_gcsi_gdsc, transcript_stability_ggdsc_ccle)
ciri_features <- rbind(ciri_stability_gcsi_ccle, ciri_stability_gcsi_gdsc, ciri_stability_ggdsc_ccle)
circ_features <- rbind(circ_stability_gcsi_ccle, circ_stability_gcsi_gdsc, circ_stability_ggdsc_ccle)
cfnd_features <- rbind(cfnd_stability_gcsi_ccle, cfnd_stability_gcsi_gdsc, cfnd_stability_ggdsc_ccle)
fcrc_features <- rbind(fcrc_stability_gcsi_ccle, fcrc_stability_gcsi_gdsc, fcrc_stability_ggdsc_ccle)


save(#gene_features, transcript_features,
     ciri_features, circ_features, cfnd_features, fcrc_features,
     file = "../results/data/temp/circ_feature_influence.RData")

############################################################
# Format results for plotting
############################################################

# function to format results for plotting
format_features <- function(features) {
  #rename datasets
  features <- features %>% 
    mutate(dataset = recode(dataset, "gcsi_ccle" = "gCSI/CCLE", "gcsi_gdsc" = "gCSI/GDSC", "gdsc_ccle" = "GDSC/CCLE"))

  #rename features
  features <- features %>% 
    mutate(feature = recode(feature, "median" = "Median Expression", "gc" = "GC%", 
                            "n_exon" = "Number of Exons", "length" = "Length"))
  features$factor <- factor(features$feature, 
                     levels=c('Median Expression', 'GC%', 'Number of Exons', 'Length'))

  # compute difference
  features$difference <- with(features, abs(baseline - permuted))

  # function from http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }

  lm_res <- data_summary(features, varname = "difference", groupnames = c("factor", "dataset"))

  return(lm_res)
}

gene_res <- format_features(gene_features)
transcript_res <- format_features(transcript_features)
ciri_res <- format_features(ciri_features)
circ_res <- format_features(circ_features)
cfnd_res <- format_features(cfnd_features)
fcrc_res <- format_features(fcrc_features)



############################################################
# Plot feature influence
############################################################

plot_features <- function(res) {

  p <- ggplot(data = res, aes(x=factor, y=difference, fill = dataset)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=difference-sd, ymax=difference+sd), width=.2, position=position_dodge(.9)) +
  ylab(expression('MSE - MSE'['permuted'])) + xlab("\nFeature") + 
  scale_fill_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("gCSI/CCLE","gCSI/GDSC","GDSC/CCLE"), values = c("#392C57", "#51C7AD", "#3670A0")) + 
  scale_x_discrete(limits = c('Median Expression', 'GC%', 'Number of Exons', 'Length'),
                   labels=c('Median\nExpression', 'GC%', 'Number\nof Exons', 'Length')) + 
  theme_classic() + scale_y_continuous(limits = c(0, 0.032), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 15), 
        legend.key.size = unit(0.6, 'cm'),
        axis.text.x = element_text(size=12, vjust = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 17), 
        legend.position = "bottom") + 
  guides(guide_legend(title.vjust = 2.8))

  return(p)
}
 

png("../results/figures/figure8/gene_features.png", width=275, height=150, units='mm', res = 600, pointsize=80)
plot_features(gene_res)
dev.off()

png("../results/figures/figure8/transcript_features.png", width=275, height=150, units='mm', res = 600, pointsize=80)
plot_features(transcript_res)
dev.off()

png("../results/figures/figure8/ciri_features.png", width=275, height=150, units='mm', res = 600, pointsize=80)
plot_features(ciri_res)
dev.off()

png("../results/figures/figure8/circ_features.png", width=275, height=150, units='mm', res = 600, pointsize=80)
plot_features(circ_res)
dev.off()

png("../results/figures/figure8/cfnd_features.png", width=275, height=150, units='mm', res = 600, pointsize=80)
plot_features(cfnd_res)
dev.off()

png("../results/figures/figure8/fcrc_features.png", width=275, height=150, units='mm', res = 600, pointsize=80)
plot_features(fcrc_res)
dev.off()