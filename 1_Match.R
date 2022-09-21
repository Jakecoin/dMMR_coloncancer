library(MatchIt)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/MergeData/")
meta <- read.csv("data/meta.csv", header=T, row.names=1)

# QC
meta <- meta[!is.na(meta$QC),]

# dmmr CRC
table(meta$Group)
table(meta$dMMR.0.否.1.是.)

# Location
table(meta$Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum.)
# Colon dMMR
table(meta$dMMR.0.否.1.是. == 1 & meta$Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum.!=2)
table(meta$dMMR.0.否.1.是. == 0 & meta$Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum.!=2)
table(meta$dMMR.0.否.1.是. %in% c(0,1) & meta$Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum.!=2)

#
meta$dmmr = NA
meta$dmmr[meta$dMMR.0.否.1.是. == 1 & meta$Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum.!=2] = "dMMR"
meta$dmmr[meta$dMMR.0.否.1.是. == 0 & meta$Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum.!=2] = "pMMR"
meta$dmmr[meta$Group == "CTRL"] = "CTRL"
table(meta$dmmr)
meta <- meta[!is.na(meta$dmmr),]
write.csv(meta, "data/meta_dmmr_full.csv")
#

meta <- read.csv("data/meta_dmmr_full.csv", header=T, row.names=1)
###
metaCRC <- meta[!(meta$dmmr %in% "CTRL"), ]
metaCRC$dmmr = factor(metaCRC$dmmr, levels = c("pMMR", "dMMR"))

out <- matchit(factor(dmmr) ~ Age + Gender, data = metaCRC, ratio = 2)
out
summary(out)

# out$match.matrix
# plot(out, type = "jitter")
x1 <- rownames(match.data(out))
# write.csv(match.data(out), "meta_matched.csv")

# CTRL
metaCTRL <- meta[meta$dmmr != "pMMR", ]
table(metaCTRL$dmmr)
dim(metaCTRL)
metaCRC$dmmr = factor(metaCRC$dmmr, levels = c("CTRL", "dMMR"))

out <- matchit(factor(dmmr) ~ Age + Gender, data = metaCTRL,
               ratio = 2)
out
summary(out)
x2 <- rownames(match.data(out))

meta <- meta[unique(c(x1, x2)),]
dim(meta)
table(meta$dmmr)

write.csv(meta, "data/meta_dmmr.csv")
# x <- read.csv("data/meta_dmmr.csv", row.names=1)
