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

countsTable <- read.table("salmon.gene.counts.f2.txt", header=TRUE, row.names = 1)
head(countsTable)
countsTableRound <-round(countsTable)
conds <- read.delim("F2.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ treatment)
dim(dds)
dds

dds <- DESeq(dds)



#list results generated
resultsNames(dds)
#[1] "Intercept"          "treatment_OA_vs_AM" "treatment_OW_vs_AM"
resF2OA_OW <- results(dds, contrast=c("treatment","OA","OW"), alpha=0.05)
resF2OA_OW <- resF2OA_OW[order(resF2OA_OW$padj),]
head(resF2OA_OW)
# log2 fold change (MLE): treatment OA vs OW 
# Wald test p-value: treatment OA vs OW 
# DataFrame with 6 rows and 6 columns
# baseMean   log2FoldChange             lfcSE             stat
# <numeric>        <numeric>         <numeric>        <numeric>
#   TRINITY_DN2010_c2_g1  2408.15398825839 2.54288220662202 0.191446357569624 13.2824789089928
# TRINITY_DN10564_c0_g1 318.632292659625 4.98252609540689 0.460600409674381 10.8174591050174
# TRINITY_DN2349_c0_g1  326.856059660863 4.64919564060459 0.445124941602375  10.444698119742
# TRINITY_DN2369_c0_g1  1942.52434666965 2.58418114828298 0.249264603032994 10.3672206837203
# TRINITY_DN7057_c0_g3  739.350568013438 2.07432642775278 0.203295231170866 10.2035173958869
# TRINITY_DN19898_c0_g1 525.126222717363 12.7846011803633  1.26187344570572 10.1314448163328
# pvalue                 padj
# <numeric>            <numeric>
#   TRINITY_DN2010_c2_g1  2.92546545271113e-40 2.39013452951952e-35
# TRINITY_DN10564_c0_g1 2.84557051124151e-27 1.16242978169471e-22
# TRINITY_DN2349_c0_g1  1.54947570128604e-25 4.21979047569236e-21
# TRINITY_DN2369_c0_g1  3.49543426448447e-25 7.13951187106613e-21
# TRINITY_DN7057_c0_g3  1.91219174417291e-24 3.12455955381342e-20
# TRINITY_DN19898_c0_g1 4.00684564243774e-24 5.45605493054676e-20
summary(resF2OA_OW)
resF2OA_OW <- resF2OA_OW[!is.na(resF2OA_OW$padj),]
degsF2OA_OW <- row.names(resF2OA_OW[resF2OA_OW$padj < 0.05,])

resF2OW <- results(dds,contrast=c("treatment","OW","AM"), alpha=0.05)
resF2OW <- resF2OW[order(resF2OW$padj),]
head(resF2OW)
# log2 fold change (MLE): treatment OW vs AM 
# Wald test p-value: treatment OW vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean    log2FoldChange             lfcSE              stat
# <numeric>         <numeric>         <numeric>         <numeric>
#   TRINITY_DN2010_c2_g1  2408.15398825839 -2.16734088056492 0.213546342211053 -10.1492765370005
# TRINITY_DN2349_c0_g1  326.856059660863 -4.80611810153229 0.489074714342944 -9.82696091330167
# TRINITY_DN10564_c0_g1 318.632292659625 -4.80721294464806 0.505692633355962 -9.50619532016044
# TRINITY_DN10337_c0_g1 208.961789109006  3.34895454034037 0.361168609949336   9.2725515121875
# TRINITY_DN28881_c0_g1 230.149399274648 -4.00683870438734 0.433271474104776 -9.24787100897019
# TRINITY_DN7057_c0_g3  739.350568013438 -2.03655954447517 0.225621790201134 -9.02643110250853
# pvalue                 padj
# <numeric>            <numeric>
#   TRINITY_DN2010_c2_g1  3.33828896285619e-24   2.503549807694e-19
# TRINITY_DN2349_c0_g1  8.61802299498734e-23 3.23154317254538e-18
# TRINITY_DN10564_c0_g1 1.97763536961932e-21 4.94375881815336e-17
# TRINITY_DN10337_c0_g1 1.81745250600778e-20 3.40749626720133e-16
# TRINITY_DN28881_c0_g1 2.29008486331284e-20 3.43489828648293e-16
# TRINITY_DN7057_c0_g3  1.77361748797232e-19 2.21687405850807e-15
summary(resF2OW)
resF2OW <- resF2OW[!is.na(resF2OW$padj),]
degsF2OW <- row.names(resF2OW[resF2OW$padj < 0.05,])

resF2OA <- results(dds, contrast=c("treatment","OA","AM"), alpha=0.05)
resF2OA <- resF2OA[order(resF2OA$padj),]
head(resF2OA)
# log2 fold change (MLE): treatment OA vs AM 
# Wald test p-value: treatment OA vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean    log2FoldChange            lfcSE              stat
# <numeric>         <numeric>        <numeric>         <numeric>
#   TRINITY_DN45760_c0_g1 44.1386056485727  21.2603398907622 2.96968375377244  7.15912590482229
# TRINITY_DN33041_c0_g1 32.8037031509363  20.7678327638551 3.02346203924104  6.86889152048634
# TRINITY_DN18365_c3_g1 46.5275742731698  19.9377006191047 2.93246867856059  6.79894750960856
# TRINITY_DN26457_c0_g1 59.3076207452932 -9.06341246312011 1.35206541163287 -6.70338312417472
# TRINITY_DN5857_c2_g1  48.2215814625406  20.5898791225979 3.05249984678778  6.74525148437442
# TRINITY_DN66575_c1_g1 32.3838269546315  20.7303395458281 3.09800398904819  6.69151480085638
# pvalue                 padj
# <numeric>            <numeric>
#   TRINITY_DN45760_c0_g1 8.11930908762156e-13 6.27038002218838e-08
# TRINITY_DN33041_c0_g1 6.47026827826389e-12 2.49842939296882e-07
# TRINITY_DN18365_c3_g1 1.05386202755922e-11 2.71292188881144e-07
# TRINITY_DN26457_c0_g1 2.03648590854349e-11 2.84291916243381e-07
# TRINITY_DN5857_c2_g1  1.52762196453197e-11 2.84291916243381e-07
# TRINITY_DN66575_c1_g1 2.20872157437754e-11 2.84291916243381e-07
summary(resF2OA)
resF2OA <- resF2OA[!is.na(resF2OA$padj),]
degsF2OA <- row.names(resF2OA[resF2OA$padj < 0.05,])

#looking at specific gene compared across treatment
plotCounts(dds, gene="TRINITY_DN45760_c0_g1", intgroup="treatment") #from F2 top: AM vs OA
plotCounts(dds, gene="TRINITY_DN2010_c2_g1", intgroup="treatment") #from F2 top: AM vs OW
plotCounts(dds, gene="TRINITY_DN2010_c2_g1", intgroup="treatment") #from F2 top: OW vs OA

plotCounts(dds, gene="TRINITY_DN29519_c0_g1", intgroup="treatment") #from F4 top: AM vs OWA
plotCounts(dds, gene="TRINITY_DN67785_c0_g1", intgroup="treatment") #from F4 top: AM vs OA
plotCounts(dds, gene="TRINITY_DN34045_c0_g1", intgroup="treatment") #from F4 top: AM vs OW

plotCounts(dds, gene="TRINITY_DN42020_c1_g1", intgroup="treatment") #form F11 top

#making venni diagram
library(eulerr)

# Total
length(degsF2OA)  # 379
length(degsF2OW)  # 1472
length(degsF2OA_OW)  #1912 

# Intersections
length(intersect(degsF2OA,degsF2OW))  # 195
length(intersect(degsF2OA,degsF2OA_OW))  # 137
length(intersect(degsF2OW,degsF2OA_OW))  # 832

intEL <- intersect(degsF2OA,degsF2OW)
length(intersect(degsF2OA_OW,intEL)) # 2
# Number unique
379-137-195-2 #45
1472-195-832-2 #443
1912-137-832-2 #941

fit1 <- euler(c("OA"=45, "OW"=443, "OA vs OW"=941, "OA&OW"=195, "OW&OA vs OW"=137, "OA&OA vs OW"=832, "OA&OW&OA vs OW"=2))

plot(fit1,  lty = 1:3, quantities = TRUE)
              
plot(fit1, quantities = TRUE, fill = "transparent",lty = 1:3,labels = list(font = 4))
              
#heatmap
library(pheatmap)

vsd <- vst(dds, blind = FALSE)

topgenes <- head(rownames(resF2OW),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(dds)[,c("treatment")])
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resF2OA),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(dds)[,c("treatment")])
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)

topgenes <- head(rownames(resF2OA_OW),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(dds)[,c("treatment")])
df <- as.data.frame(conds)
pheatmap(mat, annotation_col=df)