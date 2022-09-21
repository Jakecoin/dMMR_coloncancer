library(Maaslin2)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
# library(tidyverse)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData/")

outlier_process <- function(data){
  temp <- apply(data, 1, function(x){
    q <- quantile(x)
    iqr <- q[4]-q[2]
    median <- median(x)
    x[x > as.numeric(q[4] + 3 * iqr) | x < as.numeric(q[2] - 3 * iqr)] <- NA # 1.5 or 3
    return(x)
  })
  return(t(temp))
}

filegroup = "dmmr"
df_input_metadata <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)
dim(df_input_metadata)

target <- c("species", "met", "ko")
save_path <- c(glue("Result/{filegroup}/species_maaslin2/"), 
               glue("Result/{filegroup}/met_maaslin2/"), 
               glue("Result/{filegroup}/ko_maaslin2/"))
group.test <- c("CTRL", "pMMR", "dMMR")
# df_input_data <- outlier_process(df_input_data)

lapply(1:3, function(index){
  if(index == 1) {
    df_input_data <- read.table("data/tax_fraction_S.txt", row.names = 1)
    df_input_data <- df_input_data * 100000
  } else if (index == 2) {
    df_input_data <- read.csv("data/metabolites.csv", row.names = 1)# %>% na.omit()
  } else {
    df_input_data <- read.table("data/ko.txt", row.names = 1)
  }
  matchedL <- match_two_df(df_input_data, df_input_metadata, way = "col-row")
  df_input_data <- matchedL$df1
  df_input_metadata <- matchedL$df2
  
  # df_input_data <- outlier_process(df_input_data)
  fit_data <- Maaslin2(
    input_data = df_input_data,
    input_metadata = df_input_metadata,
    output = save_path[index],
    fixed_effects = c(filegroup, "Gender", "Age"), #"BMI",
    reference = c(filegroup, "CTRL"),
    max_significance = 0.2,
    min_prevalence = 0.1,
    min_abundance = 1)
  
  maaslin2_all_results <- fit_data$results # Save results table
  maaslin2_results <- maaslin2_all_results %>% filter(metadata == filegroup) %>% arrange(pval) # Discard covariate associations
  maaslin2_results$qval <- p.adjust(maaslin2_results$pval, method = 'BH') # FDR correction using 'BH'
  write.table(maaslin2_results, glue("{save_path[index]}{target[index]}_maaslin2_results.txt"),
              quote = FALSE)
})

lapply(1:3, function(index){
  if(index == 1) {
    df_input_data <- read.table("data/species.txt", row.names = 1)
    df_input_data <- df_input_data * 10000
    # abundance = apply(df_input_data, 1, mean)
    # table(abundance > 1)
    # df_input_data <- df_input_data[abundance > 10,]
  } else if (index == 2) {
    df_input_data <- read.csv("data/metabolites.csv", row.names = 1)# %>% na.omit()
  } else {
    df_input_data <- read.table("data/ko.txt", row.names = 1)
  }
  matchedL <- match_two_df(df_input_data, df_input_metadata, way = "col-row")
  df_input_data <- matchedL$df1
  df_input_metadata <- matchedL$df2
  idy = df_input_metadata$dmmr %in% group.test[2:3] ########
  
  fit_data <- Maaslin2(
    input_data = df_input_data[,idy],
    input_metadata = df_input_metadata[idy,],
    output = save_path[index],
    fixed_effects = c(filegroup, "Gender", "Age"),
    reference = c(filegroup, group.test[3]), #########
    max_significance = 0.2,
    min_prevalence = 0.1,
    min_abundance = 1)
  
  maaslin2_all_results <- fit_data$results # Save results table
  maaslin2_results <- maaslin2_all_results %>% filter(metadata == filegroup) %>% arrange(pval) # Discard covariate associations
  maaslin2_results$qval <- p.adjust(maaslin2_results$pval, method = 'BH') # FDR correction using 'BH'
  write.table(maaslin2_results, glue("{save_path[index]}dmmr_{target[index]}_maaslin2_results.txt"),
              quote = FALSE)
})
