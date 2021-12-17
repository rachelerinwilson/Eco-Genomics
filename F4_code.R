
setwd("~/Documents/GitHub/hudsonicaProject381/data")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(hexbin)
library(eulerr)
library(pheatmap)
countsTable <- read.table("salmon.gene.counts.matrix", header=TRUE, row.names=1)
head(countsTable)

dim(countsTable) # 130580     12

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

conds <- read.delim("ahud_samples_R.txt", header=T, stringsAsFactors = T, row.names = 1)
head(conds)

### Visualizing Data

#check how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
# just over 20 million reads per sample
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound), cex.names = 0.5, las=3,ylim=c(0,20000000)) #cex.names changes size, las changes orientation
abline(h=mean(colSums(countsTableRound)), col="blue",lwd=2)

#The average number counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #6076.078
median(rowSums(countsTableRound)) #582
#median much lower than mean; expected with distribution with long tail

apply(countsTableRound, 1, mean) #across rows
apply(countsTableRound, 2, mean) #across columns
hist(apply(countsTableRound, 1, mean), xlim=c(0,1000), breaks=10000)
hist(apply(countsTableRound, 1, mean), xlim=c(50000,150000),ylim=c(0,10), breaks=100)

##########################################################################

### create a DESeq object and define experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData=countsTableRound, 
                              colData=conds, 
                              design=~treatment+generation) #[cond]:[cond] = interaction

#+treatment:generation Error in checkFullRank(modelMatrix)
# uneven dimensions

dim(dds) #130580     38

### Filter out genes with too few reads - keep reads with average > 10 reads per sample

dds <- dds[rowSums(counts(dds)) > 380] #10 * 38 samples
dim(dds) #77724     38

### Run DESeq Model to test for differential gene expression
dds <- DESeq(dds)

resultsNames(dds)
# [1] "Intercept"            "treatment_OA_vs_AM"   "treatment_OW_vs_AM"  
# [4] "treatment_OWA_vs_AM"  "generation_F11_vs_F0" "generation_F2_vs_F0" 
# [7] "generation_F4_vs_F0"

vsd <- vst(dds, blind=F)

data <- plotPCA(vsd, intgroup=c("treatment","generation"),returnData=T)
percentVar <- round(100 *attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=treatment,shape=generation)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()
# notes: seems to be significant effect of interaction?

### TREATMENT EFFECT ###
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ generation + treatment)
dds <- DESeq(dds, test="LRT", reduced=~generation)
resultsNames(dds)
# [1] "Intercept"               
# [2] "line_combined_vs_ambient"
# [3] "environment_HH_vs_AA"   

resTreat <- results(dds, alpha = 0.05)
resTreat <- resTreat[order(resTreat$padj),] #ordering rows by adjusted p-value
head(resTreat)

summary(resTreat)
# out of 130550 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 123, 0.094%
# LFC < 0 (down)     : 109, 0.083%
# outliers [1]       : 12073, 9.2%
# low counts [2]     : 48769, 37%
# (mean count < 13)

resTreat <- resTreat[!is.na(resTreat$padj),] #filters data table to exclude NA results

degsTreat <- row.names(resTreat[resTreat$padj < 0.05,]) #differentially expressed genes from the environment
length(degsTreat) 
#232 genes differentially expressed

# TREATMENT F4 ONLY #
countsTable <- read.table("salmon.gene.counts.f4.txt", header=TRUE, row.names=1)
head(countsTable)

dim(countsTable) # 130580     12

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

conds <- read.delim("F4.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ treatment)
dim(dds)
dds

dds <- DESeq(dds)
vsd <- vst(dds, blind=F)

data <- plotPCA(vsd, intgroup=c("treatment"),returnData=T)
percentVar <- round(100 *attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Differential expression
resultsNames(dds)

#Combined vs ambient
resF4_Combo <- results(dds, alpha=0.05)
resF4_Combo <- resF4_Combo[order(resF4_Combo$padj),] #ordering rows by adjusted p-value
head(resF4_Combo)

summary(resF4_Combo)
# out of 93222 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 304, 0.33%
# LFC < 0 (down)     : 339, 0.36%
# outliers [1]       : 2724, 2.9%
# low counts [2]     : 3336, 3.6%

resF4_Combo <- resF4_Combo[!is.na(resF4_Combo$padj),]
degsCombo <- row.names(resF4_Combo[resF4_Combo$padj < 0.05,])
length(degsCombo) 
#612

#OA vs ambient
resF4_OA <- results(dds,contrast = c("treatment", "AM", "OA"), alpha=0.05)
resF4_OA <- resF4_OA[order(resF4_OA$padj),] #ordering rows by adjusted p-value
head(resF4_OA)

summary(resF4_OA)
# out of 93222 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 197, 0.21%
# LFC < 0 (down)     : 312, 0.33%
# outliers [1]       : 2724, 2.9%
# low counts [2]     : 8342, 8.9%

resF4_OA <- resF4_OA[!is.na(resF4_OA$padj),]
degsOA <- row.names(resF4_OA[resF4_OA$padj < 0.05,])
length(degsOA) 
#472

#OW vs ambient
resF4_OW <- results(dds,contrast = c("treatment", "AM", "OW"), alpha=0.05)
resF4_OW <- resF4_OW[order(resF4_OW$padj),] #ordering rows by adjusted p-value
head(resF4_OW)

summary(resF4_OW)
# out of 93222 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 185, 0.2%
# LFC < 0 (down)     : 257, 0.28%
# outliers [1]       : 2724, 2.9%
# low counts [2]     : 11635, 12%
# (mean count < 8)
resF4_OW <- resF4_OW[!is.na(resF4_OW$padj),]
degsOW <- row.names(resF4_OW[resF4_OW$padj < 0.05,])
length(degsOW) 
#404

## VENN DIAGRAM ##
# Total
length(degsOA)  # 472
length(degsOW)  # 404
length(degsCombo)  # 612

# Intersections
length(intersect(degsOA,degsOW))  # 153
length(intersect(degsOA,degsCombo))  # 144
length(intersect(degsCombo,degsOW))  # 160

intEL <- intersect(degsOA,degsOW)
length(intersect(degsCombo,intEL)) # 108

# Number unique
472-144-153-108 # 67
404-153-160-108 #-17
612-144-160-108 # 200

# unique, present in 2, present in 3
fit1 <- euler(c("OA" = 67, "OW" = 0, "Combined" = 200, 
                "OA&OW" = 153, "OA&Combined" = 144, "OW&Combined" = 160, 
                "OA&OW&Combined" = 108))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

#looking at specific gene compared across treatment
plotCounts(dds, gene="TRINITY_DN45760_c0_g1", intgroup="treatment") #from F2 top: AM vs OA
plotCounts(dds, gene="TRINITY_DN2010_c2_g1", intgroup="treatment") #from F2 top: AM vs OW
plotCounts(dds, gene="TRINITY_DN0_c0_g1", intgroup="treatment") #from F2 top: OW vs OA

plotCounts(dds, gene="TRINITY_DN29519_c0_g1", intgroup="treatment") #from F4 top: AM vs OWA
plotCounts(dds, gene="TRINITY_DN67785_c0_g1", intgroup="treatment") #from F4 top: AM vs OA
plotCounts(dds, gene="TRINITY_DN34045_c0_g1", intgroup="treatment") #from F4 top: AM vs OW

plotCounts(dds, gene="TRINITY_DN42020_c1_g1", intgroup="treatment") #form F11 top

