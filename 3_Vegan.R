library(stringr)
library(devtools)
library(ade4)
library(ggplot2)	#用于 ggplot2 作图
library(ggsignif)
library(doBy) 	#用于分组统计
library(ggalt)	#用于绘制拟合曲线
library(ggpubr)
library(dplyr)
library(scales)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(tibble)
library(vegan)	#用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(glue)

filegroup = "dmmr"
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")

mat <- read.table("data/species.txt")
groupinfo <- read.csv(glue("data/meta_{filegroup}_full.csv"), row.names = 1)
mat <- mat[,rownames(groupinfo)]

x <- round(t(mat * 10000))
alpha.index <- data.frame(richness = rowSums(x > 0), chao1 = estimateR(x)[3, ], ace = estimateR(x)[5, ],
                          shannon = diversity(x, index = 'shannon'), simpson = diversity(x, index = 'simpson'),
                          pielou = diversity(x, index = 'shannon') / log(estimateR(x)[1, ], exp(1)), 
                          goods_coverage = 1 - rowSums(x == 1) / rowSums(x),
                          group = groupinfo[[filegroup]])

table(rownames(alpha.index) == colnames(mat))
manual_color_vector = c("#276FBF", "#FF9B71", "#E84855")

gginput <- reshape2::melt(alpha.index, id.vars = c("group"), measure.vars = c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage"))
gginput$group <- factor(gginput$group, levels = c("CTRL", "pMMR", "dMMR")) ###########

unique(gginput$variable)
my_comparison <- combn(as.character(unique(groupinfo[[filegroup]])), 2, simplify=FALSE)

p <- ggplot(gginput, aes(x=group, y=value, fill=group))+
    geom_violin(trim=FALSE, aes(linetype=NA)) +
    geom_boxplot(width = 0.25, outlier.size = 0.25) +
    # geom_point(position = position_jitterdodge(),size=0.3)+
    stat_compare_means(comparisons = my_comparison,
                       method = "t.test")+
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5)) +
    scale_fill_manual(values = manual_color_vector) +
  facet_wrap(.~variable, ncol = 4, nrow = 2, scales="free_y")

p

filename <- glue("Result/{filegroup}/Alpha_index_boxplot.pdf")
pdf(filename, width = 8, height = 8)
p
dev.off()
