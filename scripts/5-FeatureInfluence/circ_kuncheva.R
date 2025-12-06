# load libraries
suppressPackageStartupMessages({
    library(genefu)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
})


############################################################
# Load in data
############################################################

load("../results/data/temp/circ_stability.RData")       # ciri: 30, circ: 66, cfnd: 84, fcrc: 503


############################################################
# Compute kuncheva index across quantiles
############################################################

# specify quantiles
quantiles <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35)

#quantiles <- c(0.15, 0.20, 0.25, 0.30, 0.35)

# function to compute kuncheva index across quantiles
compute_kuncheva <- function(stability_df, id_label) {

    # rename id_label to transcript_id
    colnames(stability_df)[which(colnames(stability_df) == id_label)] <- 'transcript_id'

    gcsi_ccle <- stability_df[,c("transcript_id","gCSI/CCLE")]
    gcsi_gdsc <- stability_df[,c("transcript_id","gCSI/GDSC2")]
    gdsc_ccle <- stability_df[,c("transcript_id","GDSC2/CCLE")]

    # initialize dataframe to store results
    kuncheva_df <- data.frame(matrix(nrow=length(quantiles), ncol = 2))
    colnames(kuncheva_df) <- c("stable","unstable")
    rownames(kuncheva_df) <- paste0("q",quantiles)

    for (q in quantiles){
    
        # compute stable overlap 
        gcsi_ccle_stable <- gcsi_ccle[order(gcsi_ccle$'gCSI/CCLE', decreasing = TRUE),][1:floor(nrow(stability_df)*q),]
        gcsi_gdsc_stable <- gcsi_gdsc[order(gcsi_gdsc$'gCSI/GDSC2', decreasing = TRUE),][1:floor(nrow(stability_df)*q),]
        gdsc_ccle_stable <- gdsc_ccle[order(gdsc_ccle$'GDSC2/CCLE', decreasing = TRUE),][1:floor(nrow(stability_df)*q),]
        
        fsets <- list()
        fsets[["gcsi_ccle"]] <- gcsi_ccle_stable$transcript_id
        fsets[["gcsi_gdsc"]] <- gcsi_gdsc_stable$transcript_id
        fsets[["gdsc_ccle"]] <- gdsc_ccle_stable$transcript_id
        
        # compute Kuncheva Index
        stable <- stab.fs(fsets=fsets, N=nrow(stability_df), method="kuncheva")
        kuncheva_df[ paste0("q",q), "stable"] <- stable

        # compute unstable overlap 
        gcsi_ccle_unstable <- gcsi_ccle[order(gcsi_ccle$'gCSI/CCLE', decreasing = FALSE),][1:floor(nrow(stability_df)*q),]
        gcsi_gdsc_unstable <- gcsi_gdsc[order(gcsi_gdsc$'gCSI/GDSC2', decreasing = FALSE),][1:floor(nrow(stability_df)*q),]
        gdsc_ccle_unstable <- gdsc_ccle[order(gdsc_ccle$'GDSC2/CCLE', decreasing = FALSE),][1:floor(nrow(stability_df)*q),]
        
        fsets <- list()
        fsets[["gcsi_ccle"]] <- gcsi_ccle_unstable$transcript_id
        fsets[["gcsi_gdsc"]] <- gcsi_gdsc_unstable$transcript_id
        fsets[["gdsc_ccle"]] <- gdsc_ccle_unstable$transcript_id
        
        # compute Kuncheva Index
        unstable <-  stab.fs(fsets=fsets, N=nrow(stability_df), method="kuncheva")
        kuncheva_df[ paste0("q",q), "unstable"] <- unstable
    }

    # formating for plotting
    kuncheva_df$quantile <- as.character(quantiles)
    kuncheva <- melt(kuncheva_df)

    return(kuncheva)
}

ciri_kuncheva <- compute_kuncheva(ciri_stability, "circ_id")
circ_kuncheva <- compute_kuncheva(circ_stability, "circ_id")
cfnd_kuncheva <- compute_kuncheva(cfnd_stability, "circ_id")
fcrc_kuncheva <- compute_kuncheva(fcrc_stability, "circ_id")



############################################################
# Plot distribution of kuncheva index
############################################################
  
# function to plot distribution of Kuncheva index
plot_kuncheva <- function(melted_df, title) {

    p <- ggplot(melted_df, aes(x=quantile, y=value, group=variable), size = 3) +
    geom_line(aes(color=variable), size = 1)+
    geom_point(aes(color=variable)) + xlab("Quantile") + ylab("Kuncheva Index") + labs(color='Stability') + 
    scale_color_manual(values=c("#8B7B96", "#71A2B6")) + 
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
            text = element_text(size = 17), 
            legend.key.size = unit(0.5, 'cm'),
            axis.text=element_text(size=15), 
            axis.text.x = element_text(size=15), 
            axis.text.y = element_text(size=15), 
            legend.title = element_text(size=15), 
            legend.text = element_text(size=15)) +
    ggtitle(title)

    return(p)
}

# plot line graph
png("../results/figures/figure8/kuncheva.png", width=250, height=200, units='mm', res = 600, pointsize=80)
ggarrange(plot_kuncheva(ciri_kuncheva, "CIRI2"),                # max q: 0.2
          plot_kuncheva(circ_kuncheva, "CIRCexplorer2"),        # max q: 0.1
          plot_kuncheva(cfnd_kuncheva, "circRNA_finder"),       # max q: 0.15
          plot_kuncheva(fcrc_kuncheva, "find_circ"),            # max q: 0.1
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()
