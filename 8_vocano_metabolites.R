library(stringr)
library(devtools)
library(ggplot2)	#用于 ggplot2 作图
library(ggsignif)
library(scales)
library(ggsci)
library(ggrepel)
library(data.table)
library(glue)

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

rawmat = read.table(glue("{filename}_VIP.txt"), row.names = 1)
rawmat = rawmat[rownames(mat),]

p.threshold = 0.05
log2fc.threshold = 0.5

diff_output <- data.frame(
  name = rownames(mat),
  abundance = apply(mat, 1, median),
  VIP = rawmat,
  metabolite = metname[rownames(mat),]$MS2Metabolite,
  pvalue = apply(mat, 1, function(x){wilcox.test(unlist(x[meta[[filegroup]] == "dMMR"]), unlist(x[meta[[filegroup]] == "pMMR"]))$p.value}),
  log2fc = apply(mat,1,function(x){log2(median(x[meta[[filegroup]] == "dMMR"] + 1)/(median(x[meta[[filegroup]] == "pMMR"] + 1)))})
)  %>% arrange(pvalue) %>% mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
  mutate(label=ifelse(log2fc > log2fc.threshold, "UP", ifelse(log2fc < (-log2fc.threshold), "DOWN", "Nosig"))) %>%
  mutate(TPplotlabel = ifelse((pvalue < p.threshold & abs(log2fc) > log2fc.threshold) & VIP > 1, metabolite, NA)) %>%
  mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel), 0, 1),sep = "")))

write.csv(diff_output, glue("{filename}.p{p.threshold}.log2fc{log2fc.threshold}.csv"))

x.axis.lim <- round(max(abs(diff_output$log2fc))+1)
y.axis.lim <- (max(-log10(diff_output$pvalue)))

pdf(glue("{filename}.volcano.p{p.threshold}.log2fc{log2fc.threshold}.pdf"), width = 15, height = 15)
ggplot(data=as.data.frame(diff_output), aes(x=log2fc, y=-log10(pvalue))) +
  geom_point(aes(color = TPplotcolor), alpha = 0.8) + # size = log10(concentration)
  scale_color_manual(labels = c("Nosig0"="No Significance","UP1"="Significantly Up","UP0"="Up","DOWN1"="Significantly Down","DOWN0"="Down"),
                     values = c("UP0"='#FAC0AE',"DOWN0"='#9BCFF0',"UP1"='#FA2311',"DOWN1"='#6175DB',"Nosig0" ='gray')) +
  labs(x="log2 Fold Change",  y="-log10 pvalue") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.position = c(0.12, 0.75)) +
  theme(legend.title = element_blank(), 
        legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent')) +
  geom_hline(yintercept = -log10(p.threshold), linetype=4, color = 'black', size = 0.5) +
  geom_vline(xintercept = c(-log2fc.threshold, log2fc.threshold), linetype = 4, color = 'black', size = 0.5) +
  geom_text_repel(data=diff_output[!duplicated(diff_output$metabolite), ], 
                  aes(x = log2fc, y = -log10(pvalue), label = TPplotlabel), size=4, colour = "black",
                  max.overlaps = 20) +
  xlim(-x.axis.lim, x.axis.lim) + ##控制横坐标长度
  ylim(1, y.axis.lim)+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
dev.off()
