library(stringr)
library(devtools)
library(ggbiplot)
library(picante)	#用于计算 PD_whole_tree，若不计算它就无需加载。
library(ggplot2)	#用于 ggplot2 作图
library(ggsignif)
library(doBy) 	#用于分组统计
library(ggalt)	#用于绘制拟合曲线
library(ggpubr)
library(dplyr)
library(scales)
library(ggsci)
library(data.table)
library(readxl)
library(glue)
library(ImageGP)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")

filter.spe <- function(x) {
  return(x[grep(x, pattern = "CAG|sp\\.|unclassified|_bacterium", invert = TRUE)])
}

myggbox <- function(mat, selected) {
  # selected = downfeature
  s.mat <- mat
  s.mat <- s.mat[rownames(s.mat) %in% selected,]
  s.abundance = apply(s.mat, 1, median)
  rownames(s.mat)[order(s.abundance, decreasing = T)]
  
  if(index == 2) {
    temp = metname[rownames(s.mat), 1]
    !duplicated(temp)
    s.mat = s.mat[!duplicated(temp), ]
    rownames(s.mat) = metname[rownames(s.mat), 1]
  }
  
  ggplot_input <- as.data.frame(reshape::melt(as.matrix(s.mat)) %>% dplyr::select(item=1,sample=2,expr=3) %>% 
                                  mutate(group=meta[sample,filegroup], expr=as.numeric(expr)))
  
  ggplot_input$group <- factor(ggplot_input$group, levels = group.test, labels = group.test)
  ggplot_input$item <- factor(ggplot_input$item, levels = rownames(s.mat)[order(s.abundance, decreasing = T)])
  
  library(ggplot2)
  col <- c("#276FBF", "#FF9B71", "#E84855")
  p <- ggplot(ggplot_input, aes(x=item, y=log2(expr + 1), fill=group))+
    geom_boxplot(width=0.8, outlier.size = 0.5) +
    scale_fill_manual(values = col) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    # geom_point(position = position_jitterdodge(), size=0.1)+
    stat_compare_means(aes(group = group),
                       label = "p.signif",
                       method = "wilcox.test",
                       hide.ns = T)+
    theme(
      legend.position = "right",
      panel.background = element_blank(),
      panel.border = element_rect(color = 'black', fill = "transparent"),
      axis.text = element_text(color = 'black',size = 6, 
                               face = 'plain'),
      axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    ) +
    labs(x="", y="log2(Relative Abundance)")
  return(p)
}

# filter.f <- function(dat){
#   Num <- ncol(dat)
#   SD <- apply(dat,1,sd)
#   num_0 <- apply(dat, 1, function(x) length(which(x == 0)))
#   ave_abun <- apply(dat,1,mean)
#   tmp <- cbind(dat,data.frame(SD),data.frame(num_0),data.frame(ave_abun))
#   colnames(tmp)[(ncol(dat)+1):(ncol(dat)+3)] <- c("sd","count0",'aveabun')
#   #dat_filter <- tmp %>% filter(count0 <= as.numeric(Num*0.9) & sd >0) 
#   dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0) & (tmp$aveabun > 0.00001),]
#   #dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0),]
#   return(dat_filter[,1:ncol(dat)])
# }

filegroup = "dmmr"
group.test <- c("CTRL", "pMMR", "dMMR")
meta <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)

index = 1
target = c("species", "metabolites", "kogene")
filename = glue("Result/{filegroup}/boxplot_{target[index]}")

i = "S"
if(index == 1) {
  # mat <- read.table("data/species.txt", row.names = 1)
  mat <- read.table(glue("data/tax_{i}.txt"))
  # mat <- mat * 100000
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

group = meta[[filegroup]]

# abundance threshold
abundance.threshold = 10
abundance = apply(mat, 1, median)
table(abundance >= abundance.threshold)
mat <- mat[abundance >= abundance.threshold,]

diff_output <- data.frame(
  name = rownames(mat), # metname[rownames(mat),1]
  abundance = apply(mat, 1, mean),
  pvalue = apply(mat, 1, function(x){wilcox.test(unlist(x[group == group.test[3]]), unlist(x[group == group.test[2]]))$p.value}),
  log2fc = apply(mat,1,function(x){log2(median(x[group == group.test[3]] + 1)/(median(x[group == group.test[2]] + 1)))})
)  %>% arrange(pvalue) %>% mutate(fdr=p.adjust(pvalue, method = "BH"))

p.threshold = 0.05
log2fc.threshold = 0.5
diff_output <- diff_output[order(diff_output$abundance, decreasing = T),]
write.table(diff_output, glue("{filename}_{i}_diff.txt"))

upfeature <- diff_output[diff_output$pvalue < p.threshold & diff_output$log2fc > 0,]$name %>% na.omit()
downfeature <- diff_output[diff_output$pvalue < p.threshold & diff_output$log2fc <= 0,]$name %>% na.omit()
upfeature
downfeature

write.table(upfeature, glue("{filename}_{i}_up.txt"))
write.table(downfeature, glue("{filename}_{i}_down.txt"))

p1 <- myggbox(mat, upfeature)
p2 <- myggbox(mat, downfeature)

pdf(glue("{filename}_{p.threshold}_{i}_up.pdf"), onefile = T, width = length(upfeature)/2, height = 6) #length(selected)/2
p1
dev.off()

pdf(glue("{filename}_{p.threshold}_{i}_down.pdf"), onefile = T, width = length(downfeature)/2, height = 6) #length(selected)/2
p2
dev.off()

# ggarrange(p1, p2, align = "hv")
