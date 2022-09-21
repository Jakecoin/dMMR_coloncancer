# rm(list=ls())
library(WGCNA)
library(igraph)
library(psych)
library(impute)
library(corrplot)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData")

filegroup = "dmmr"
group.test <- c("CTRL", "pMMR", "dMMR")
meta <- read.csv(glue("data/meta_{filegroup}.csv"), row.names = 1)

i = "S"
otu <- read.table(glue("data/tax_{i}.txt"))
ev <- read.csv("data/metabolites.csv", row.names = 1) %>% na.omit()
metname <- ev[,c(558:568)]
ko <- read.table("data/ko.txt", row.names = 1)

matchedL <- match_two_df(otu, meta, way = "col-row")
meta <- matchedL$df2
matchedL <- match_two_df(ev, meta, way = "col-row")
meta <- matchedL$df2
matchedL <- match_two_df(ko, meta, way = "col-row")
meta <- matchedL$df2

table(meta[[filegroup]])

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

evup <- metname[evup, 1]
evdown <- metname[evdown, 1]

otusig = data.frame(feature = c(otuup, otudown),
                    group = c(rep("Enriched", length(otuup)),
                              rep("Depleted", length(otudown)))
                    )
evsig = data.frame(feature = c(evup, evdown),
                    group = c(rep("Enriched", length(evup)),
                              rep("Depleted", length(evdown)))
)
evsig <- evsig[!duplicated(evsig$feature),]
kosig = data.frame(feature = c(koup, kodown),
                   group = c(rep("Enriched", length(koup)),
                             rep("Depleted", length(kodown)))
)

otusig$group2 <- "species"
evsig$group2 <- "metabolites"
kosig$group2 <- "kogene"
# kosig$group2 <- kosig$subtype # type
# kosig <- subset(kosig, select = c("feature", "group", "group2"))

feature.temp <- rbind(otusig, evsig) %>% rbind(kosig)
rownames(feature.temp) <- feature.temp$feature

temp <- rbind(otu, ev)
temp <- rbind(temp, ko)

# start
rvalue <- 0.6
ko.rvalue <- 0.6
pvalue <- 0.005
point.del <- 1
withname <- TRUE
pdf(glue("Result/{filegroup}/correlation_r{rvalue}_rko{ko.rvalue}_p{pvalue}_del{point.del}_{withname}.pdf"),
    width = 24, height = 8)

par(mfrow=c(1, 3))

layouts = c("layout_with_fr", "layout_with_kk", "layout_with_dh",
            "layout_with_gem", "layout_as_star", "layout_as_tree",
            "layout_in_circle", "layout_on_grid")

lapply(1:3, function(i){
  OTU <- t(temp)
  idy <- meta[[filegroup]] %in% group.test[i]
  table(idy)
  OTU <- OTU[idy,]

  occor = corAndPvalue(OTU)
  
  occor.r = occor$cor # 取相关性矩阵R值
  occor.p = occor$p # 取相关性矩阵p值
  
  table(occor.p < pvalue & abs(occor.r) > rvalue)
  occor.r[occor.p >= pvalue | abs(occor.r) <= rvalue] = 0 
  
  idx1 <- 1:nrow(otu)
  idx2 <- (nrow(otu)+1):(nrow(otu)+nrow(ev)) #ev
  idx3 <- (nrow(otu)+nrow(ev)+1):nrow(temp) #ko
  idx4 <- (nrow(otu)+1):nrow(temp) #ev+ko
  # occor.r[idx1, idx1] = 0
  occor.r[idx2, idx2] = 0
  occor.r[idx3, idx3] = 0
  occor.r[idx4, idx4] = 0
  #################################################
  for(k in idx1) {
    for(j in idx3) {
      occor.r[k, j] = ifelse(abs(occor.r[k, j]) < ko.rvalue, 0, occor.r[k, j])
      occor.r[j, k] = ifelse(abs(occor.r[j, k]) < ko.rvalue, 0, occor.r[j, k])
    }
  }

  igraph = graph_from_adjacency_matrix(occor.r, mode = "undirected", weighted=TRUE, diag=FALSE)

  if (point.del > 0) {
    bad.vs = V(igraph)[igraph::degree(igraph) < point.del]
    igraph = igraph::delete.vertices(igraph, bad.vs)
  }
  igraph.weight = E(igraph)$weight
  E(igraph)$weight = NA
  set.seed(123)
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#D04539",ifelse(E.color<0, "#3B94BA", "grey"))
  E(igraph)$color = as.character(E.color)
  V(igraph)$size = 5
  
  # set vertices color
  rownames(feature.temp) <- feature.temp$feature
  igraph.col = feature.temp[V(igraph)$name,]
  
  colist <- c("#A31621", "#FCF7F8", "#AB9F9D", "#4E8098", "#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6", "#086788", "#5B5941", "#E9D2F4", "#340068", "#D88C9A")
  if(length(colist) <= length(unique(igraph.col$group2))) {
    colist <- rep(colist, length.out = length(unique(igraph.col$group2)))
  } else {
    colist <- colist[1:length(unique(igraph.col$group2))]
  }
  
  names(colist) <- unique(igraph.col$group2)
  ver.col <- factor(igraph.col$group2, levels = names(colist), labels = colist)
  V(igraph)$color = as.character(ver.col)
  set.seed(123)
  
  l=do.call(layouts[4], list(igraph))
  if(withname) {
    plot(igraph, main=paste(group.test[i], " network", sep = ""),
       vertex.frame.color="black", edge.arrow.size = 1, # vertex.label = NA,
       edge.lty=1, edge.curved=FALSE, vertex.label.dist = 0)
  } else {
    plot(igraph, main=paste(group.test[i], " network", sep = ""),
         vertex.frame.color="black", vertex.label = NA, edge.arrow.size = 1,
         edge.lty=1, edge.curved=FALSE)
  }
  # legend("right", legend = names(colist), fill = colist)
  # x <- as_data_frame(igraph, what="edges")
  # y <- as_data_frame(igraph, what="vertices")
})
dev.off()

#####################################################################################################################################
rvalue <- 0.3
ko.rvalue <- 0.6
pvalue <- 0.005
point.del <- 1
withname <- FALSE

i = 3

pdf(glue("Result/{filegroup}/correlation_{group.test[i]}_r{rvalue}_rko{ko.rvalue}_p{pvalue}_del{point.del}_{withname}.pdf"),
    width = 8, height = 8)

layouts = c("layout_with_fr", "layout_with_kk", "layout_with_dh",
            "layout_with_gem", "layout_as_star", "layout_as_tree",
            "layout_in_circle", "layout_on_grid")

if(i == 3){
  feature.TEMP <- feature.temp[c(otuup, evup, koup),]
  idy <- meta[[filegroup]] %in% group.test[3]
  OTU <- t(temp)[idy,c(otuup, evup, koup)]
} else if (i == 2) {
  feature.TEMP <- feature.temp[c(otudown, evdown, kodown),]
  idy <- meta[[filegroup]] %in% group.test[2]
  OTU <- t(temp)[idy, c(otudown, evdown, kodown)]
}
{
  occor = corAndPvalue(OTU)
  
  occor.r = occor$cor # 取相关性矩阵R值
  occor.p = occor$p # 取相关性矩阵p值
  
  table(occor.p < pvalue & abs(occor.r) > rvalue)
  occor.r[occor.p > pvalue | abs(occor.r) < rvalue] = 0 
  
  # idx1 <- 1:nrow(otu)
  # idx2 <- (nrow(otu)+1):(nrow(otu)+nrow(ev)) #ev
  # idx3 <- (nrow(otu)+nrow(ev)+1):nrow(temp) #ko
  # idx4 <- (nrow(otu)+1):nrow(temp) #ev+ko
  # # occor.r[idx1, idx1] = 0
  # occor.r[idx2, idx2] = 0
  # occor.r[idx3, idx3] = 0
  # occor.r[idx4, idx4] = 0
  #################################################
  # for(k in idx1) {
  #   for(j in idx3) {
  #     occor.r[k, j] = ifelse(abs(occor.r[k, j]) < ko.rvalue, 0, occor.r[k, j])
  #     occor.r[j, k] = ifelse(abs(occor.r[j, k]) < ko.rvalue, 0, occor.r[j, k])
  #   }
  # }
  
  igraph = graph_from_adjacency_matrix(occor.r, mode = "undirected", weighted=TRUE, diag=FALSE)
  
  if (point.del > 0) {
    bad.vs = V(igraph)[igraph::degree(igraph) < point.del]
    igraph = igraph::delete.vertices(igraph, bad.vs)
  }
  igraph.weight = E(igraph)$weight
  E(igraph)$weight = NA
  set.seed(123)
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#D04539",ifelse(E.color<0, "#3B94BA", "grey"))
  E(igraph)$color = as.character(E.color)
  V(igraph)$size = 5
  
  # set vertices color
  igraph.col = feature.TEMP[V(igraph)$name,]
  
  colist <- c("#A31621", "#FCF7F8", "#AB9F9D", "#4E8098", "#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6", "#086788", "#5B5941", "#E9D2F4", "#340068", "#D88C9A")
  if(length(colist) <= length(unique(igraph.col$group2))) {
    colist <- rep(colist, length.out = length(unique(igraph.col$group2)))
  } else {
    colist <- colist[1:length(unique(igraph.col$group2))]
  }
  names(colist) <- unique(igraph.col$group2)
  ver.col <- factor(igraph.col$group2, levels = names(colist), labels = colist)
  V(igraph)$color = as.character(ver.col)
  set.seed(123)
  
  l=do.call(layouts[4], list(igraph))
  if(withname) {
    plot(igraph, main=paste(group.test[i], " network", sep = ""),
         vertex.frame.color="black", edge.arrow.size = 1, # vertex.label = NA,
         edge.lty=1, edge.curved=FALSE, vertex.label.dist = 0)
  } else {
    plot(igraph, main=paste(group.test[i], " network", sep = ""),
         vertex.frame.color="black", vertex.label = NA, edge.arrow.size = 1,
         edge.lty=1, edge.curved=FALSE)
  }
  # legend("right", legend = names(colist), fill = colist)
}
dev.off()

