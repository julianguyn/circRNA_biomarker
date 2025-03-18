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
stability_pair <- function(stability) {
  #todo
}

#Make a copy for each dataset pair
transcript_stability_gcsi_ccle <- transcript_stability
transcript_stability_gdsc_ccle <- transcript_stability
transcript_stability_gcsi_gdsc <- transcript_stability


#Isolate MOST stable/unstable transcripts (15% quantile)
#gCSI/CCLE
stable <- transcript_stability_gcsi_ccle[order(-transcript_stability_gcsi_ccle$gcsi_ccle_spearman),]
stable <- stable[1:9037,] #top 15%
unstable <- transcript_stability_gcsi_ccle[order(transcript_stability_gcsi_ccle$gcsi_ccle_spearman),]
unstable <- unstable[1:9037,]
transcript_stability_gcsi_ccle <- rbind(stable,unstable)

#gCSI/GDSC
stable <- transcript_stability_gcsi_gdsc[order(-transcript_stability_gcsi_gdsc$gcsi_gdsc_spearman),]
stable <- stable[1:9037,]
unstable <- transcript_stability_gcsi_gdsc[order(transcript_stability_gcsi_gdsc$gcsi_gdsc_spearman),]
unstable <- unstable[1:9037,]
transcript_stability_gcsi_gdsc <- rbind(stable,unstable)

#GDSC/CCLE
stable <- transcript_stability_gdsc_ccle[order(-transcript_stability_gdsc_ccle$gdsc_ccle_spearman),]
stable <- stable[1:9037,]
unstable <- transcript_stability_gdsc_ccle[order(transcript_stability_gdsc_ccle$gdsc_ccle_spearman),]
unstable <- unstable[1:9037,]
transcript_stability_gdsc_ccle <- rbind(stable,unstable)


#Extract features + stability from df
transcript_stability_gcsi_ccle <- transcript_stability_gcsi_ccle[,c('gdsc_median','gc','n_exon','single_read_mappability','multi_read_mappability','length','three_prime_bias','five_prime_bias','gcsi_ccle_spearman')]
transcript_stability_gcsi_gdsc <- transcript_stability_gcsi_gdsc[,c('ccle_median','gc','n_exon','single_read_mappability','multi_read_mappability','length','three_prime_bias','five_prime_bias','gcsi_gdsc_spearman')]
transcript_stability_gdsc_ccle <- transcript_stability_gdsc_ccle[,c('gcsi_median','gc','n_exon','single_read_mappability','multi_read_mappability','length','three_prime_bias','five_prime_bias','gdsc_ccle_spearman')]

#shuffle data
transcript_stability_gcsi_ccle <- transcript_stability_gcsi_ccle[sample(1:nrow(transcript_stability_gcsi_ccle)),]
transcript_stability_gcsi_gdsc <- transcript_stability_gcsi_gdsc[sample(1:nrow(transcript_stability_gcsi_gdsc)),]
transcript_stability_gdsc_ccle <- transcript_stability_gdsc_ccle[sample(1:nrow(transcript_stability_gdsc_ccle)),]


########################################################
########### Permutation Feature Importance #############
########################################################


################# 
# Reduce size of data.frame for Code Ocean, comment out following 3 lines for full result
transcript_stability_gcsi_ccle <- transcript_stability_gcsi_ccle[1:500,]
transcript_stability_gcsi_gdsc <- transcript_stability_gcsi_gdsc[1:500,]
transcript_stability_gdsc_ccle <- transcript_stability_gdsc_ccle[1:500,]
#################


## Run Linear Regression
runLinearReg <- function(x, label) {
  
  # Create dataframe to store feature scores
  feature_score <- data.frame(matrix(ncol=4, nrow = 0))
  colnames(feature_score) <- c("feature", "baseline", "permuted", "dataset")
  
  # Change colnames of inputted dataframe
  colnames(x) <- c("median", colnames(x)[2:8], "spearman")
  
  #10-fold CV
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 5) #change to 10 for full results
  
  
  # Loop through each feature
  for (i in 1:length(colnames(x)[1:8])) {
    
    feature <- colnames(x)[i]
    
    # get baseline score
    gbmFit <- train(spearman ~ ., data = x, 
                    method = "lm", 
                    trControl = fitControl,
                    verbose = FALSE)
    #Save MSE score
    baseline_score <- (gbmFit$results$RMSE)^2
    
    # Permute each feature 20,000 times
    for (k in 1:2) { # change to 20,000 for full results
      
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


lm_gcsi_ccle <- runLinearReg(transcript_stability_gcsi_ccle, "gcsi_ccle")
lm_gcsi_gdsc <- runLinearReg(transcript_stability_gcsi_gdsc, "gcsi_gdsc")
lm_gdsc_ccle <- runLinearReg(transcript_stability_gdsc_ccle, "gdsc_ccle")

lm_feature_results <- rbind(lm_gcsi_ccle, lm_gcsi_gdsc, lm_gdsc_ccle)



###############
#### PLOTS ####
###############

load("C:/Users/julia/Desktop/BHK Lab/permutations/h4h/h4h_lm_results.RData")

#rename datasets
lm_feature_results <- lm_feature_results %>% 
  mutate(dataset = recode(dataset, "gcsi_ccle" = "gCSI/CCLE", "gcsi_gdsc" = "gCSI/GDSC", "gdsc_ccle" = "GDSC/CCLE"))

#rename features
lm_feature_results <- lm_feature_results %>% 
  mutate(feature = recode(feature, "median" = "Median Expression", "gc" = "GC%", "n_exon" = "Number of Exons", "length" = "Length",
                          "single_read_mappability" = "Single Read Mappability", "multi_read_mappability" = "Multi Read Mappability",
                          "three_prime_bias" = "Three Prime Bias", "five_prime_bias" = "Five Prime Bias"))


lm_feature_results$factor <- factor(lm_feature_results$feature, 
                                    levels=c('Median Expression', 'GC%', 'Number of Exons', 'Length', 'Single Read Mappability', 'Multi Read Mappability', 'Three Prime Bias', 'Five Prime Bias'))


#### Feature Influence ####

lm_feature_results$difference <- with(lm_feature_results, abs(baseline - permuted))

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
lm_res <- data_summary(lm_feature_results, varname = "difference", groupnames = c("factor", "dataset"))

p1 <- ggplot(data = lm_res, aes(x=factor, y=difference, fill = dataset)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=difference-sd, ymax=difference+sd), width=.2, position=position_dodge(.9)) +
  ylab(expression('MSE - MSE'['permuted'])) + xlab("\nFeature") + 
  scale_fill_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("gCSI/CCLE","gCSI/GDSC","GDSC/CCLE"), values = c("#392C57", "#51C7AD", "#3670A0")) + 
  scale_x_discrete(limits = c('Median Expression', 'GC%', 'Number of Exons', 'Length', 'Single Read Mappability', 'Multi Read Mappability', 'Three Prime Bias', 'Five Prime Bias'),
                   labels=c('Median\nExpression', 'GC%', 'Number\nof Exons', 'Length', 'Single Read\nMappability','Multi Read\nMappability', 'Three Prime\nBias','Five Prime\nBias')) + 
  theme_classic() + scale_y_continuous(limits = c(0, 0.032), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 15), 
        legend.key.size = unit(0.6, 'cm'),
        axis.text.x = element_text(size=12, vjust = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 17), 
        legend.position = "bottom") + 
  guides(guide_legend(title.vjust = 2.8)) 

png("../results/figure3.png", width=275, height=150, units='mm', res = 600, pointsize=80)
p1
dev.off()