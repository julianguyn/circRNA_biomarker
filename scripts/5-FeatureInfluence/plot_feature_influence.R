# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
})

set.seed(200)

source("utils/palettes.R")

############################################################
# Load in elastic net spearman results
############################################################

# get files
path <- "../results/data/feature_influence/"
files <- list.files(path, recursive = TRUE, pattern = ".*en.csv")

# initiate dataframe for storing results
df <- data.frame(matrix(nrow=0, ncol=0))

# get results
for (file in files) {
  res <- read.csv(paste0(path, "/", file))
  res$dataset <- sub("/.*", "", file)
  res$pair <- sub(".*/", "", sub("/en.csv", "", file))
  df <- rbind(df, res)
}

# format results
df$fold <- factor(df$fold, levels = c(1:5))
df$pair <- factor(
  df$pair,
  levels = c("gcsi_ccle", "gcsi_gdsc", "gdsc_ccle"),
  labels = c("gCSI & CCLE", "gCSI & GDSC2", "CCLE & GDSC2")
)
df$dataset <- factor(df$dataset, levels = names(dataset_pal))
df$cat <- paste(df$dataset, df$pair, sep = "-")


############################################################
# Summarize spearman across folds
############################################################

# get average and max
res <- df %>%
  group_by(cat) %>%
  summarise(
    avg_spearman = mean(spearman, na.rm = TRUE),
    ci_spearman  = sd(spearman, na.rm = TRUE)/ sqrt(5) * 1.96,
    max_spearman = spearman[which.max(abs(spearman))],
    .groups = "drop"
  )

# format results again
res$label <- gsub("-.*", "", res$cat)
res$pair <- gsub(".*-", "", res$cat)
res$label <- factor(res$label, levels = names(dataset_pal))
res$pair <- factor(res$pair, levels = names(prop_pal[2:4]))

res <- res[order(res$label, res$pair),] |> as.data.frame()

############################################################
# Plot model results (average across folds)
############################################################

# set fill colours for CI past zero
toPlot <- res %>%
  mutate(
    crosses_zero = (avg_spearman - ci_spearman <= 0 & avg_spearman + ci_spearman >= 0),
    fill = ifelse(crosses_zero, "Non-Significant", as.character(pair))
  )
toPlot$fill <- factor(toPlot$fill, levels = c(names(prop_pal[2:4]), "Non-Significant"))


png("../results/figures/suppfig4/en_spearman.png", width=10, height=5, units='in', res = 600, pointsize=80)
ggplot(toPlot, aes(x = pair, y = avg_spearman, fill = fill)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg_spearman-ci_spearman, ymax=avg_spearman+ci_spearman), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c(prop_pal, "gray")) +
  facet_wrap(.~label, nrow = 1) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.size = unit(0.5, 'cm'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
  ) +
  labs(fill = "Dataset Pair", x = "Dataset Pair", y = "Average Spearman") 
dev.off()