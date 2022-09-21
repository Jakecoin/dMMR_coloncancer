library(corrplot)
library(dplyr)
require(readxl)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(tibble)

library(igraph)
library(psych)
library(impute)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")

filegroup = "dmmr"
group.test <- c("CTRL", "pMMR", "dMMR")
meta <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)

index = 1
target = c("species", "metabolites", "kogene")
filename = glue("Result/{filegroup}/boxplot_{target[index]}")

i = "S"
otu <- read.table(glue("data/tax_{i}.txt"))
ev <- read.csv("data/metabolites.csv", row.names = 1) %>% na.omit()
metname <- mat[,c(558:568)]
ko <- read.table("data/ko.txt", row.names = 1)

matchedL <- match_two_df(otu, meta, way = "col-row")
meta <- matchedL$df2
matchedL <- match_two_df(ev, meta, way = "col-row")
meta <- matchedL$df2
matchedL <- match_two_df(ko, meta, way = "col-row")
meta <- matchedL$df2

group = meta[[filegroup]]

otuup <- read.table(glue("Result/{filegroup}/alteration/boxplot_species_S_up.txt"))[,1] # species_heatmap_sig_0518.txt
evup <- read.table(glue("Result/{filegroup}/alteration/boxplot_metabolites_up.txt"))[,1]
koup <- read.table(glue("Result/{filegroup}/alteration/boxplot_kogene_up.txt"))[,1]

otudown <- read.table(glue("Result/{filegroup}/alteration/boxplot_species_S_down.txt"))[,1] # species_heatmap_sig_0518.txt
evdown <- read.table(glue("Result/{filegroup}/alteration/boxplot_metabolites_down.txt"))[,1]
kodown <- read.table(glue("Result/{filegroup}/alteration/boxplot_kogene_down.txt"))[,1]

otu <- otu[c(otuup, otudown), rownames(meta)]
ev <- ev[c(evup, evdown), rownames(meta)]
ko <- ko[c(koup, kodown), rownames(meta)]

temp = metname[rownames(ev), 1]
ev = ev[!duplicated(temp), ]
rownames(ev) = metname[rownames(ev), 1]

cormat <- corr.test(t(otu), t(ev), method="spearman", use="complete.obs")
# write.csv(cormat, file="spearman_correlation_.csv", quote=F, row.names=T)
# mattest <- cor.mtest(cormat, conf.level = 0.95)
# namesmat1 = c(otuup[1], otuup[length(otuup)], metname[evup[1], 1], metname[evup[length(evup)], 1])
# namesmat2 = c(otudown[1], otudown[length(otudown)], metname[evdown[1], 1], metname[evdown[length(evdown)], 1])
# namesmat = rbind(namesmat1, namesmat2)

pdf(glue("Result/{filegroup}/OTUvsEV_spearman_correlation.pdf"), width = nrow(ev)/6, height = nrow(otu)/6)
corrplot(-cormat$r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
         p.mat = cormat$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5, insig = 'label_sig',
         col = COL2('RdYlBu')) # %>% corrRect(namesMat = namesmat)
dev.off()

cormat <- corr.test(t(otu), t(ko), method="spearman", use="complete.obs")
pdf(glue("Result/{filegroup}/OTUvsKO_spearman_correlation.pdf"), width = nrow(ko)/6, height = nrow(otu)/6)
corrplot(-cormat$r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
         p.mat = cormat$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5, insig = 'label_sig',
         col = COL2('RdYlBu'))
dev.off()

cormat <- corr.test(t(ev), t(ko), method="spearman", use="complete.obs")
pdf(glue("Result/{filegroup}/EVvsKO_spearman_correlation.pdf"), width = nrow(ko)/6, height = nrow(ev)/6)
corrplot(-cormat$r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
         p.mat = cormat$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5, insig = 'label_sig',
         col = COL2('RdYlBu'))
dev.off()

# library(reshape2)
# nodes <- as.data.frame(mat) %>% dplyr::mutate(item = rownames(.), abun = rowMeans(.)) %>%
#                  dplyr::select(item, abun)
# 
# edges <- cbind(cormat %>% reshape2::melt(), (mattest$p %>% reshape2::melt())$value)
# colnames(edges) <- c("from", "to", "value", "pvalue")
# edges <- edges[edges$pvalue<0.05 & edges$value!=1 & abs(edges$value) > 0.4,]
# 
# e.pvalue <- edges$pvalue
# edges <- edges[,-4]
# 
# ig <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed = FALSE)
# tg <- tidygraph::as_tbl_graph(ig) %>% 
#   tidygraph::activate(nodes) %>% 
#   dplyr::mutate(label=name)
# 
# # set seed
# set.seed(12345)
# 
# tg %>%
#   ggraph(layout = "fr") +
#   geom_edge_arc(colour= "gray50",
#                 lineend = "round",
#                 strength = .1,
#                 alpha = .1) +
#   geom_node_text(aes(label = name), 
#                  repel = TRUE, 
#                  point.padding = unit(0.2, "lines"), 
#                  colour="gray10") +
#   theme_graph(background = "white") +
#   guides(edge_width = FALSE,
#          edge_alpha = FALSE)
# 
# v.size <- log(nodes$abun)+2
# # inspect
# v.size
# 
# e.weight <- ifelse(abs(edges$value)>0.3, edges$value, NA)
# e.color <- ifelse(e.weight > 0, "Positive", "Negative")
# 
# species <- rownames(mat)[1:22]
# metabolite <- rownames(mat)[23:28]
# Family <- c(rep("Species", 22), rep("metabolite", 6))
# 
# col_fun <- colorRamp2(c(-0.8,0,0.8), c("#0F8B8D", "white", "#A8201A"))
# 
# pdf("network.pdf")
# 
# # edge size shows frequency of co-occurrence
# tg %>%
#   ggraph(layout = "fr") +
#   geom_edge_arc(#colour= "gray50",
#                 lineend = "round",
#                 strength = .1,
#                 aes(edge_width = e.weight,
#                     edge_alpha = e.weight,
#                     colour = e.color)) +
#   geom_node_point(size=sqrt(v.size)*2, 
#                   aes(colour=Family, shape=Family)) +
#   geom_node_text(aes(label = name), 
#                  repel = TRUE, 
#                  point.padding = unit(0.2, "lines"), 
#                  size=4, 
#                  colour="gray10") +
#   scale_edge_width(range = c(0, 2.5)) +
#   scale_edge_alpha(range = c(0, .3)) +
#   # scale_edge_color_manual(
#   #   limits = as.factor(layout$name),
#   #   values = cols_f(nrow(layout))
#   # ) +
#   theme_graph(background = "white") +
#   theme(legend.position = "top") 
#   # guides(edge_width = FALSE,
#   #        edge_alpha = FALSE)
# 
# dev.off()
# 
# ######################################
# # FC_col_fun = colorRamp2(c(-7, 0, 7), c("#3BB273", "white", "#7768AE"))
# 
# # ha = rowAnnotation(log2FC = anno_simple(metmat_diff, col = FC_col_fun), 
# #                    show_annotation_name = FALSE
# #                    #annotation_name_gp = gpar(fontsize = 6)
# # )
# col = colorRamp2(c(-0.6, 0, 0.6), c("#18578B", "white", "#8B1818"))
# 
# pdf("Ob_stricted_Spearman_correlation.pdf",width = 8,height=5)
# 
# ht <- Heatmap(htmat, name="Spearman correlation", col = col,
#               column_names_side = "top", column_names_gp = gpar(fontsize = 8),
#               column_names_rot = 45,
#               row_names_gp = gpar(fontsize = 8),
#               # row_km = 4, column_km = 4,
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.rect(x = x, y = y, width = width, height = height,
#                           gp = gpar(col = "#E0E0E0", fill = NA))
#                 grid.circle(x = x, y = y, r = abs(htmat[i, j])/2 * min(unit.c(width, height)),
#                             gp = gpar(fill = col(htmat[i, j]), col=NA))
#               },
#               layer_fun = function(j, i, x, y, width, height, fill) {
#                 v = pindex(htmat, i, j)
#                 l = abs(v) > 0.4
#                 grid.text(sprintf("%.2f", v[l]), x[l], y[l], gp = gpar(fontsize = 5))
#               }
# )
# ht
# # log2FC = Legend(title = "log2FC", col_fun = FC_col_fun, at = c(-7, 0, 7), labels = c("-7", "0", "7"))
# 
# # draw(ht, annotation_legend_list = list(log2FC))
# 
# dev.off()
# 
# # 
# # idmax <- apply(matrix,1,function(x){max(abs(x))})>0.2
# # View(matrix[idmax,])
# # dim(matrix[idmax,])
# # 
# # cormat <- matrix
# # 
# # mat.test <- cor.mtest(cormat, conf.level = 0.95)
# # # 
# # pdf("Young34.MvsN.spearman_correlation.pdf", width=20, height=20)
# # corrplot(matrix, order = "hclust", addrect = FALSE,
# #          #method = colour,
# #          col = colorRampPalette(c("#C4315D", "white", "#60B2E5"))(100),
# #          is.corr = FALSE, col.lim = c(-0.35,0.35),
# #          #p.mat = mat.test$p, insig = "label_sig", sig.level = c(0.001,0.01,0.05), pch.cex=0.9, pch.col="black", 
# #          )
# # dev.off()
# 
# # p.mat = mat.test$p, sig.level = 0.05, insig = "label_sig", sig.level = c(0.001,0.01,0.05), pch.cex=0.9, pch.col="black", 
# # cl.lim = c(-0.4, 0.4)
# 
# 
# #~~~~~~~~
# 
# #O
# matrix <- cor(t(met.sp.matrix.O), t(sp.matrix.O),method="spearman",use="complete.obs")
# 
# dim(matrix)
# matrix[1:5,1:5]
# write.csv(matrix,file="O.MvsN.spearman_correlation.csv",quote=F,row.names=T)
# 
# metmat_diff <- met_O_diff_output[met_O_diff_output$pvalue<0.05 & abs(met_O_diff_output$log2FC)>2,]$log2FC
# 
# col = colorRamp2(c(-0.4, 0, 0.4), c("#18578B", "white", "#8B1818"))
# 
# FC_col_fun = colorRamp2(c(-6, 0, 6), c("#3BB273", "white", "#7768AE"))
# 
# ha = rowAnnotation(log2FC = anno_simple(metmat_diff, col = FC_col_fun), 
#                    show_annotation_name = FALSE
#                    #annotation_name_gp = gpar(fontsize = 6)
#                    )
# 
# pdf("test.Old25_MvsN_Spearman_correlation.pdf",width = 12,height=8)
# ht <- Heatmap(matrix, name="Spearman correlation", col = col,
#               right_annotation = ha,
#               rect_gp = gpar(type = "none"),
#               column_names_side = "top", column_names_gp = gpar(fontsize = 8),
#               column_names_rot = 45,
#               row_names_gp = gpar(fontsize = 8),
#               #row_km = 4, column_km = 4,
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.rect(x = x, y = y, width = width, height = height, 
#                           gp = gpar(col = "#E0E0E0", fill = NA))
#                 grid.circle(x = x, y = y, r = abs(matrix[i, j])/2 * min(unit.c(width, height)) * 3, 
#                             gp = gpar(fill = col(matrix[i, j]), col=NA))
#               }
#               # layer_fun = function(j, i, x, y, width, height, fill) {
#               #   v = pindex(matrix, i, j)
#               #   l = abs(v) > 0.1
#               #   grid.text(sprintf("%.1f", v[l]), x[l], y[l], gp = gpar(fontsize = 10))
#               # }
# )
# log2FC = Legend(title = "log2FC", col_fun = FC_col_fun, at = c(-6, 0, 6), labels = c("-6", "0", "6"),)
# draw(ht, annotation_legend_list = list(log2FC), heatmap_legend_side = "right", annotation_legend_side = "right")
# 
# dev.off()

