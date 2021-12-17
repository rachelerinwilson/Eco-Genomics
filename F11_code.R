#F2 data
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn) 

setwd("~/EcoGenom/Project DATA")

countsTable <- read.table("salmon.gene.counts.f11.txt", header=TRUE, row.names = 1)
head(countsTable)
countsTableRound <-round(countsTable)
conds <- read.delim("F11.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ treatment)
dim(dds)
dds

dds <- DESeq(dds)



#list results generated
resultsNames(dds)
#[1] "Intercept"           "treatment_OWA_vs_AM"
resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)
# log2 fold change (MLE): treatment OWA vs AM 
# Wald test p-value: treatment OWA vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean   log2FoldChange             lfcSE              stat
# <numeric>        <numeric>         <numeric>         <numeric>
#   TRINITY_DN42020_c1_g1 526.660949406619  3.9263567978919 0.473579501152697  8.29080817124709
# TRINITY_DN39749_c0_g1 170.923700832766 -11.020432083344  1.37257642433908 -8.02901163674762
# TRINITY_DN13962_c0_g2 164.793250240051 10.7444640669775  1.37005256196118  7.84237361783926
# TRINITY_DN30572_c0_g1 157.082399312488 10.6758243135466  1.35582453763257  7.87404565799337
# TRINITY_DN32275_c0_g1 81.4271128442866 9.72789458705983  1.31325685766298  7.40745767310991
# TRINITY_DN17545_c0_g1  111.68910352498 10.1833943070657  1.44912124183449  7.02728937585248
# pvalue                 padj
# <numeric>            <numeric>
#   TRINITY_DN42020_c1_g1 1.12481564941091e-16  9.6294342930419e-12
# TRINITY_DN39749_c0_g1 9.82610386911643e-16 4.20601463065594e-11
# TRINITY_DN13962_c0_g2 4.42108608822753e-15 9.46211897317676e-11
# TRINITY_DN30572_c0_g1 3.43353547396586e-15 9.46211897317676e-11
# TRINITY_DN32275_c0_g1 1.28743554895525e-13 2.20432139821021e-09
# TRINITY_DN17545_c0_g1 2.10584342075958e-12 2.00310166008674e-08
summary(resInt)



#looking at specific gene compared across treatment
plotCounts(dds, gene="TRINITY_DN45760_c0_g1", intgroup="treatment") #from F2 top: AM vs OA
plotCounts(dds, gene="TRINITY_DN2010_c2_g1", intgroup="treatment") #from F2 top: AM vs OW
plotCounts(dds, gene="TRINITY_DN2010_c2_g1", intgroup="treatment") #from F2 top: OW vs OA

plotCounts(dds, gene="TRINITY_DN29519_c0_g1", intgroup="treatment") #from F4 top: AM vs OWA
plotCounts(dds, gene="TRINITY_DN67785_c0_g1", intgroup="treatment") #from F4 top: AM vs OA
plotCounts(dds, gene="TRINITY_DN34045_c0_g1", intgroup="treatment") #from F4 top: AM vs OW

plotCounts(dds, gene="TRINITY_DN42020_c1_g1", intgroup="treatment") #form F11 top

resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])
length(degsInt) #1898

#making venni diagram
# library(eulerr)
# fit1 <- euler(c("AM"=, "OWA"=))
# plot(fit1,  lty = 1:3, quantities = TRUE)

#heatmap
library(pheatmap)

vsd <- vst(dds, blind = FALSE)

topgenes <- head(rownames(resInt),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(dds)[,c("treatment")])
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)
