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

two_feature_heatmap <- function(mat, upfeature, downfeature){
  s.mat <- apply(mat[c(upfeature, downfeature),], 1, scale)
  s.mat <- t(s.mat)
  colnames(s.mat) <- colnames(mat)
  
  column_ha = HeatmapAnnotation(Status = meta$dmmr,
                                col = list(Status = c("CTRL" = "#276FBF", "pMMR" = "#FF9B71", "dMMR" = "#E84855"))
  )
  
  signature = c(rep("Enriched", length(upfeature)), rep("Depleted", length(downfeature)))
  row_ha = rowAnnotation(Alteration = signature,
                         col = list(Alteration = c("Enriched" = "#FA2311", "Depleted" = "#6175DB"))
  )
  
  lim = round(max(abs(s.mat))/2)-1
  
  if(index == 2)
    rowlabels = metname[rownames(s.mat), 1]
  
  p <- Heatmap(s.mat[c(upfeature, downfeature),], name = "Z-score", border = TRUE, #rownames(s.mat) %in% 
               show_column_names = FALSE, show_row_names = TRUE,
               row_labels = rowlabels,
               top_annotation = column_ha,
               left_annotation = row_ha,
               column_split = factor(meta[[filegroup]], levels = group.test, labels = group.test),
               cluster_column_slices = FALSE,
               column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
               row_names_gp = gpar(fontsize = 10),
               row_order = c(upfeature, downfeature),
               column_names_rot = 0, column_names_centered = TRUE,
               col = colorRamp2(c(-lim,0,lim), c("#6175DB", "white", "#FA2311"))
  )
  pdf(glue("{filename}_heatmap.pdf"), width = dim(s.mat)[2]/8, height = dim(s.mat)[1]/5) # 
  print(p)
  dev.off()
  return(0)
}
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")
filegroup = "dmmr"
group.test <- c("CTRL", "pMMR", "dMMR")
meta <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)

index = 2
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
} else {
  mat <- read.table("data/ko.txt", row.names = 1)
  
}

matchedL <- match_two_df(mat, meta, way = "col-row")
mat <- matchedL$df1
meta <- matchedL$df2
if(index == 2){
  mat <- log2(mat+1)
}

filename = glue("Result/{filegroup}/{target[index]}/boxplot_{target[index]}")
upfeature = read.table(glue("{filename}_up.txt"))[,1]
downfeature = read.table(glue("{filename}_down.txt"))[,1]

# heatmap
two_feature_heatmap <- function(mat, upfeature, downfeature, subfeature = NULL){
  if(!is.null(subfeature)) {
    s.mat <- mat[rownames(metname[metname[, 2] == subfeature,]),]
    upfeature <- upfeature[upfeature %in% rownames(metname[metname[, 2] == subfeature,])]
    downfeature <- downfeature[downfeature %in% rownames(metname[metname[, 2] == subfeature,])]
  } else {
    s.mat <- mat
  }
  s.mat <- s.mat[c(upfeature, downfeature),]
  s.mat <- apply(mat[c(upfeature, downfeature),], 1, scale)
  s.mat <- t(s.mat)
  colnames(s.mat) <- colnames(mat)
  
  column_ha = HeatmapAnnotation(Status = meta[[filegroup]],
                                col = list(Status = c("CTRL" = "#276FBF", "pMMR" = "#FF9B71", "dMMR" = "#E84855"))
  )
  
  signature = c(rep("Enriched", length(upfeature)), rep("Depleted", length(downfeature)))
  row_ha = rowAnnotation(Alteration = signature,
                         col = list(Alteration = c("Enriched" = "#FA2311", "Depleted" = "#6175DB"))
  )
  
  lim = round(max(abs(s.mat))/2)-1
  
  if(index == 2)
    rowlabels = metname[rownames(s.mat), 1]
  
  p <- Heatmap(s.mat[c(upfeature, downfeature),], name = "Z-score", border = TRUE, #rownames(s.mat) %in% 
          show_column_names = FALSE, show_row_names = TRUE,
          row_labels = rowlabels,
          top_annotation = column_ha,
          left_annotation = row_ha,
          column_split = factor(meta[[filegroup]], levels = group.test, labels = group.test),
          cluster_column_slices = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 10),
          row_order = c(upfeature, downfeature),
          column_names_rot = 0, column_names_centered = TRUE,
          col = colorRamp2(c(-lim,0,lim), c("#6175DB", "white", "#FA2311"))
  )
  filepathname = ifelse(is.null(subfeature), glue("{filename}_heatmap.pdf"),
                        glue("{filename}_heatmap_{subfeature}.pdf"))
  pdf(filepathname, width = dim(s.mat)[2]/8, height = dim(s.mat)[1]/6 + 1.5) # 
  print(p)
  dev.off()
  return(0)
}


two_feature_heatmap(mat, upfeature, downfeature, subfeature = NULL) # , subfeature = "Lipids and lipid-like molecules"

for(i in unique(metname[c(upfeature, downfeature), 2])){
  if(table(metname[c(upfeature, downfeature), 2] == i)[2] > 1)
  two_feature_heatmap(mat, upfeature, downfeature, subfeature = i)
}

