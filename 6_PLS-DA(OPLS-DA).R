# PLS-DA
# https://www.blog4xiang.world/posts/3f4f7ea1.html

library(ade4)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
require(readxl)
library(ropls)

filename <- "plsda0917"

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")
filegroup = "dmmr"
group.test <- c("CTRL", "pMMR", "dMMR")
meta <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)

index = 2
target = c("species", "metabolites", "kogene")
filename = glue("Result/{filegroup}/plsda_{target[index]}")
mat <- read.csv("data/metabolites.csv", row.names = 1) %>% na.omit()
metname <- mat[,c(558:568)]

matchedL <- match_two_df(mat, meta, way = "col-row")
mat <- matchedL$df1
meta <- matchedL$df2
mat <- log2(mat+1)

idy = meta$dmmr!="CTRL"

s.mat <- as.data.frame(apply(mat, 1, function(x){scale(x)}))
rownames(s.mat) <- colnames(mat)
s.mat[1:5,1:5]

met.plsda <- opls(s.mat[idy,], meta[[filegroup]][idy], permI = 100, predI = 1, orthoI = 1) # multigroup不能设定好level , 

sample.score <- met.plsda@scoreMN %>% 
                      as.data.frame() %>%
                      mutate(o1 = met.plsda@orthoScoreMN[,1],
                              group = meta[[filegroup]][idy])#, o1 = met.plsda@orthoScoreMN[,1])
    
ggplot(sample.score, aes(p1, o1, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  labs(title = paste0("p=",met.plsda@summaryDF[8]),
       x = paste('t1(', round(met.plsda@modelDF$R2X[1]*100), '%)', sep = ''),
       y = paste('o1(', round(met.plsda@modelDF$R2X[2]*100), '%)', sep = '')) + ## 横纵坐标名字可能需要修改 round(met.plsda@modelDF$R2X[1]*100)
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) +
  scale_color_manual(values = c('#84DCC6','#EB8884')) +
  theme_bw() +
  theme(
        legend.position = "none", # right
        legend.title = element_blank(),
        legend.text = element_text(color = 'black',size = 15, face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, face = 'plain'),
        axis.title = element_text(color = 'black',size = 15, face = 'plain'),
        axis.ticks = element_line(color = 'black')
  )
    
ggsave(glue("{filename}.pdf"), width=5, height=5)
out_met <- met.plsda@vipVn[order(met.plsda@vipVn, decreasing = T)]

write.table(out_met, glue("{filename}_VIP.txt"),
            quote = FALSE, col.names = FALSE)
