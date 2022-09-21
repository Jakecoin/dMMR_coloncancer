library(lattice)
library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(tibble)
library(RColorBrewer)

library(caret)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(devtools)

library(biomaRt)
library(ComplexHeatmap)
library(circlize)

# heatmap
ko_feature_heatmap <- function(mat, feature, featuretype){
  koname <- feature
  koname <- kotext[kotext$ko %in% koname,]
  table(koname$type)
  ko.duplicated <- duplicated(koname$ko)
  koname <- koname[!ko.duplicated,]
  write.table(koname, glue("{filename}_pathway_annotation.txt"))
  
  s.mat <- mat[koname$ko,]
  s.mat <- t(apply(s.mat, 1, scale))
  colnames(s.mat) <- colnames(mat)
  
  koname$type <- factor(koname$type, levels = unique(koname$type))
  table(koname$type)
  colist <- c("#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6", "#086788", "#E9D2F4", "#5B5941", "#340068", "#D88C9A", "#000000")
  length(colist)
  length(unique(koname$type))
  colist <- rep(colist, length.out = length(unique(koname$type)))
  names(colist) <- unique(koname$type)
  
  ha <- rowAnnotation(pathway = koname$type, 
                      col = list(pathway = colist)
  )
  
  split = data.frame(type = koname$type)
  row_title_rot = 0
  
  column_ha = HeatmapAnnotation(Status = meta[[filegroup]],
                                col = list(Status = c("CTRL" = "#276FBF", "pMMR" = "#FF9B71", "dMMR" = "#E84855"))
  )
  
  lim = 1 #round(max(abs(s.mat))/2) - 1
  pdf(glue("{filename}_heatmap_pathway_{featuretype}.pdf"), width = dim(s.mat)[2]/8, height = dim(s.mat)[1]/6 + 1.5)
  print(
    Heatmap(s.mat, name = "Z-score", border = TRUE,
            left_annotation = ha,
            row_split = split, 
            top_annotation = column_ha,
            column_split = factor(meta[[filegroup]], levels = group.test, labels = group.test),
            cluster_column_slices = FALSE,
            row_title = unique(koname$type),
            row_title_rot = 0,
            row_title_gp = gpar(fontsize = 8),
            row_names_gp = gpar(fontsize = 8),
            show_column_names = FALSE,# show_row_names = FALSE,
            column_names_side = "top", 
            show_column_dend = FALSE, 
            show_row_dend = FALSE,
            col = colorRamp2(c(-lim, 0, lim), c("#6175DB", "white", "#FA2311")),
            #, row_order = rownames(s.mat)
    )
  )
  dev.off()
  return(0)
}

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")
filegroup = "dmmr"
group.test <- c("CTRL", "pMMR", "dMMR")
meta <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)

index = 3
target = c("species", "metabolites", "kogene")
filename = glue("Result/{filegroup}/boxplot_{target[index]}")

if(index == 1) {
  mat <- read.table("data/species.txt", row.names = 1)
  mat <- mat * 10000
  # abundance = apply(mat, 1, mean)
  # table(abundance > 1)
  # mat <- mat[abundance > 10,]
} else if (index == 2) {
  mat <- read.csv("data/metabolites.csv", row.names = 1) %>% na.omit()
  metname <- mat[,c(558:568)]
  mat <- log2(mat+1)
} else {
  mat <- read.table("data/ko.txt", row.names = 1)
  mat <- log2(mat+1)
}

matchedL <- match_two_df(mat, meta, way = "col-row")
mat <- matchedL$df1
meta <- matchedL$df2

filename = glue("Result/{filegroup}/{target[index]}/boxplot_{target[index]}")
upfeature = read.table(glue("{filename}_up.txt"))[,1]
downfeature = read.table(glue("{filename}_down.txt"))[,1]

kotext <- read.table("data/ko_database.txt")
colnames(kotext) <- kotext[1,]
kotext <- kotext[-1,]
kotext <- kotext[,-1]

two_feature_heatmap(mat, upfeature, downfeature, subfeature = NULL) # , subfeature = "Lipids and lipid-like molecules"

ko_feature_heatmap(mat, upfeature, "up")
ko_feature_heatmap(mat, downfeature, "down")
