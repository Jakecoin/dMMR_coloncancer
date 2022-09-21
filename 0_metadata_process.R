library(lattice);library(Formula);library(ggplot2);
library(readxl)
library(tidyverse)
library(ggpubr)
library(grid)
library(vcd)
library(stringr)
library(jstable)
library(tableone)
library(stringr)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")

meta <- read_excel("sample_map.xlsx", sheet = "Sheet1", col_names = TRUE) %>% as.data.frame()
rownames(meta) <- meta$SampleID
table(meta$Group5)

meta <- read_excel("Full_meta.xlsx", sheet = "CANCER", col_names = TRUE) %>% as.data.frame()
dim(meta)
meta = meta[,c(10,5,6,11:68)]
rownames(meta) <- meta[,1]

###
colnames(meta)[1:4] <- c("SampleID", "Gender", "Age", "Group")

meta <- within(meta, {
  Gender[Gender==0] <- "Female"
  Gender[Gender==1] <- "Male"
  Gender[Gender=="女"] <- "Female"
  Gender[Gender=="男"] <- "Male"
})

write.csv(meta, "data/meta.csv")

### Species # fraction_
meta <- read.csv("data/meta.csv", row.names = 1)
for(i in c("C", "D", "F", "G", "O", "P", "S", "S1")){
  mat <- read_excel(glue("Taxonomy_Abundance/Merge_abunance_stat_{i}.xlsx"), # S # G
                    sheet = "Sheet1", col_names = TRUE) %>% as.data.frame()
  rownames(mat) <- gsub(mat$name, pattern = "[[:punct:]]", replacement = '.')
  # tax <- mat[,1:2]
  mat <- mat[,rownames(meta)]
  write.table(mat, glue("data/tax_{i}.txt"))
}

temp = read.table(glue("data/tax_S1.txt"))
# 
for(i in c( "G", "S")){ #"C", "D", "F", "O", "P", 
  mat <- read.table(glue("data/tax_{i}.txt"))
  temp <- rbind(temp, mat)
}
write.table(temp, glue("data/tax_GS.txt"))

### Metabolites
mat <- fread("data/metabolites_raw.csv")
dim(mat)
mat[1:5,1:5]
table(mat$IsMS2ident == 1)
mat <- mat[mat$IsMS2ident == 1, ]
dim(mat)
mat$ID <- gsub(mat$ID, pattern = "[[:punct:]]", replacement = '.')
table(is.na(mat))
colnames(mat)
table(is.na(mat[1:547]))
write.csv(mat, "data/metabolites.csv", row.names = F)
# mat <- read.csv("data/metabolites.csv", row.names = 1)

### KO Entry
mat <- read_excel("Function/KEGG/4_KOEntry/KEGG_KOEntry_abund.xlsx", 
                  sheet = "Sheet1", col_names = TRUE) %>% as.data.frame()
rownames(mat) <- mat[,1]
mat <- mat[,-1]
# mat <- mat[,rownames(meta)]
write.table(mat, "data/ko.txt")
# x <- read.table("data/ko.txt")

# Input
# if(index == 1) {
#   mat <- read.table("data/species.txt", row.names = 1)
#   mat <- mat * 10000
#   # abundance = apply(mat, 1, mean)
#   # table(abundance > 1)
#   # mat <- mat[abundance > 10,]
# } else if (index == 2) {
#   mat <- read.csv("data/metabolites.csv", row.names = 1)# %>% na.omit()
# } else {
#   mat <- read.table("data/ko.txt", row.names = 1)
# }
# matchedL <- match_two_df(mat, meta, way = "col-row")
# mat <- matchedL$df1
# meta <- matchedL$df2
