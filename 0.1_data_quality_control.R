# 0_data_quality_control
library(glue)
library(stats)
library(ggplot2)
library(reshape2)
library(dplyr)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")

taxfilter <- function(mat, lowquar, highquar){
  # Bacteria
  iqr = quantile(mat[2,], c(lowquar, highquar))
  iqr
  bac.filtered <- colnames(mat)[mat[2,] < iqr[,1]]
  bac.filtered
  
  # Eukaryota
  iqr = quantile(mat[1,], c(lowquar, highquar))
  iqr
  eur.filtered <- colnames(mat)[mat[1,] > iqr[,2]]
  eur.filtered
  
  Virus
  iqr = quantile(mat[3,], c(lowquar, highquar))
  iqr
  vir.filtered <- colnames(mat)[mat[3,] > iqr[,2]]
  vir.filtered
  # # Archaea
  # iqr = quantile(mat[4,], c(lowquar, highquar))
  # iqr
  # arc.filtered <- colnames(mat)[mat[4,] > iqr[,2]]
  # arc.filtered
  
  # unique(c(bac.filtered, eur.filtered, vir.filtered))
  tax.filtered = unique(c(bac.filtered, eur.filtered, vir.filtered)) #, vir.filtered, arc.filtered))
  return(tax.filtered)
}

metagenome_name = c("D", "P", "C", "O", "F", "G", "S", "S1")

i = 1
mat = read.table(glue("data/tax_fraction_{metagenome_name[i]}.txt"))

# tmat <- t(mat)
# logmat <- t(log10(mat * 10000 + 1))
# boxplot(t(mat[,]))
# boxplot(logmat[,])

mat <- mat[,order(mat[2,], decreasing = T)]
pcm <- reshape2::melt(t(mat), id = colnames(mat)) %>% dplyr::select(sample=1, tax=2, value=3)
manual_scale_colours = c("#247BA0",  "#FF6978", "#ED9B40", "#AA8F66")

p <- ggplot(pcm, aes(x = sample, fill = tax, y = value)) +
  geom_col() +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"), 
    legend.text = element_text(colour = "black", size = 12, face = "bold"), 
    legend.title = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Relative Abundance (%)", fill = "tax") +
  scale_fill_manual(name = "Phylum", values = manual_scale_colours)

p
ggsave(plot = p, "Result/raw_stackedbarplot.pdf", width = ncol(mat)/30, height = 3)

avg.abundance = apply(mat, 1, median)
avg.abundance

summary.abundance = apply(mat, 1, summary)
summary.abundance

tax.filtered = taxfilter(mat, 0.1, 0.9)

pcm <- reshape2::melt(t(mat[, !(colnames(mat) %in% tax.filtered)]), id = colnames(mat)) %>% dplyr::select(sample=1, tax=2, value=3)
manual_scale_colours = c("#247BA0",  "#FF6978", "#ED9B40", "#AA8F66")

p <- ggplot(pcm, aes(x = sample, fill = tax, y = value)) +
  geom_col() +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"), 
    legend.text = element_text(colour = "black", size = 12, face = "bold"), 
    legend.title = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Relative Abundance (%)", fill = "tax") +
  scale_fill_manual(name = "Phylum", values = manual_scale_colours)
p
ggsave(plot = p, "Result/filtered_stackedbarplot.pdf", 
       width = ncol(mat[, !(colnames(mat) %in% tax.filtered)])/30, height = 3)

# write.table(mat[, !(colnames(mat) %in% tax.filtered)], "data/tax_fraction_{metagenome_name[i]}_filtered.txt")
selected = colnames(mat[, !(colnames(mat) %in% tax.filtered)])
meta <- read.csv(glue("data/meta.csv"), row.names = 1)

meta$QC = NA
meta[selected,]$QC = TRUE
# write.csv(meta, "data/meta.csv")
meta = meta[meta$QC == 1,]
table(meta$dMMR.0.否.1.是.)

# rnacol = c("yCRC_134", "yCRC_85", "yCRC_57", "yCRC_19", "oCRC_46", "oCRC_67")
rnacol = c("CRC67", "CRC139", "CRC165", "CRC258", "CRC270", "CRC294")
apply(mat[, colnames(mat) %in% rnacol], 1, summary)
table(selected %in% rnacol)

p <- ggplot(pcm[pcm$sample %in% rnacol,], aes(x = sample, fill = tax, y = value)) +
  geom_col() +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"), 
    legend.text = element_text(colour = "black", size = 12, face = "bold"), 
    legend.title = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Relative Abundance (%)", fill = "tax") +
  scale_fill_manual(name = "Phylum", values = manual_scale_colours)
p
