#sample trim: warming (HA) F4 = HA_F4_Rep

#navigate to directory:
$cd /data/project_data/RNAseq/rawdata/
  
  #wildcard = * says to grab anything that matches
  #fastqc FILENAME*.fq.gz --outdir=/data/project_data/RNAseq/fastqc/: putting into fastqc
  $ fastqc HA_F4* --outdir=/data/project_data/RNAseq/fastqc/
    
    #next view html file by move to local by file transfer and open by double click
    in fastqc $pwd /data/project_data/RNAseq/fastqc/ 
    $ll
  # will have to move my file from fastqc to myresults so can access 
  
  #make html file to view all FastQC for all samples be in fastqc directory. for raw data file
  $ multiqc .
  #try get into files Ihave access to
  multiqc ~/myresults
  #script for Trimmomatic
  #!/bin/bash   
  
  # Below is a script to loop through the files in the /rawdata directory, identify matches,  
  # and clean the fastq files, and direct output to /cleandata 
  
  cd /data/project_data/RNAseq/rawdata
  
  for f1 in *_1.fq.gz  
  
  do 
  
  f2=${f1%%_1.fq.gz}"_2.fq.gz"  
  
  java -classpath /data/popgen/Trimmomatic-0.39/trimmomatic-0.39.jar org.usadellab.trimmomatic.TrimmomaticPE \
  -threads 10 \
  -phred33 \
  "$f1" \
  "$f2" \
  /data/project_data/RNAseq/cleandata/"$f1"_left_clean_paired.fq \
  /data/project_data/RNAseq/cleandata/"$f1"_left_clean_unpaired.fq \
  /data/project_data/RNAseq/cleandata/"$f2"_right_clean_paired.fq \
  /data/project_data/RNAseq/cleandata/"$f2"_right_clean_unpaired.fq \
  ILLUMINACLIP:/data/popgen/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:20 \
  TRAILING:20 \
  SLIDINGWINDOW:6:20 \
  MINLEN:36 
  >> log.txt
  
  #end of script loop
  
  #ll html files
  $ll *.html
  # give count of number html files in directory
  $ll *.html | wc -1
  
  # run FastQC, get right directory, then run FastQC
  $ cd /data/project_data/RNAseq/cleandata
  $ fastqc FILENAME*.fq.gz --outdir=/data/project_data/RNAseq/fastqc/clean
  # my file HA_F4
  $ fastqc HA_F4*.fq.gz --outdir=/data/project_data/RNAseq/fastqc/clean
  # compile all all html togther in one html be where clean data is so in above directory. for clean files
  $ pwd /data/project_directory/RNAseq/fastqc/clean
  $ multiqc .
    
  #run FastQC on both before and after cleaning (Trimmomatic)

#Transcriprome Assembly: done for us though Trinity: basic command
  Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G

#Salmon to quantify transcript abundance 
  #index ref transcriptome
  cd /data/project_data/RNAseq/assembly/
    conda activate salmon
  
  salmon index -t Bridger.fasta -i hudsonica_index -p 8
  
  #in screen run salmon
  #!/bin/bash
  ######
  #
  # quantify each sample with salmon
  #
  #######
  
  # -i points to the index files already created
  # -l A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.)
  # -p 8 says uses 8 threads
  # -o indicates the directory and name of output
  # seqbias corrects for random hexamer priming
  # gcbias corrects for gcbias, but only when present.
  
  conda activate salmon
  
  for i in $(ls /data/project_data/RNAseq/cleandata | grep '.fq.gz' | cut -f 1-3 -d "_"| uniq);
  do
  
  echo "starting sample ${i}"
  #starting with only name of rep. need to pull out files
  
  read1=$(ls /data/project_data/RNAseq/cleandata | grep ${i} | grep '_1.qc.fq.gz')
  read2=$(ls /data/project_data/RNAseq/cleandata | grep ${i} | grep '_2.qc.fq.gz')
  
  salmon quant -i /data/project_data/RNAseq/assembly/hudsonica_index \
  -l A \
  -1 /data/project_data/RNAseq/cleandata/${read1} \
  -2 /data/project_data/RNAseq/cleandata/${read2} \
  -p 8  \
  --softclip \
  --seqBias \
  --gcBias \
  -o /data/project_data/RNAseq/salmon/transcripts_quant/${i}
  
  echo "sample ${i} done"
  
  done
#####Now move from command line to working in R for rest:
  
# set working directory
setwd("~/EcoGenom")

# get working dir
getwd()

## Import or install the libraries that we're likely to need
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  

#import counts matrix
countsTableF1 <- read.table("DE_counts_F1.txt", header=TRUE, row.names = 1)
countsTableF3 <- read.table("DE_counts_F3.txt", header=TRUE, row.names = 1)
dim(countsTableF1)
#[1] 24362    16
dim(countsTableF3)
#[1] 25279    16

#because DESeq2 not want decimals, salmon outputs with decimals
countsTableRoundF1 <-round(countsTableF1) 
countsTableRoundF3 <-round(countsTableF3)

#import sample discription table
condsF1 <- read.delim("RT_tonsa_F1_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
condsF3 <- read.delim("RT_tonsa_F3_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)

#STATs: see reads for each sample
colSums(countsTableRoundF1)
mean(colSums(countsTableRoundF1))
barplot(colSums(countsTableRoundF1), names.arg = colnames(countsTableRoundF1), cex.names = 0.5, las=3, ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRoundF1)), col="blue", lwd=2)

colSums(countsTableRoundF3)
mean(colSums(countsTableRoundF3))
barplot(colSums(countsTableRoundF3), names.arg = colnames(countsTableRoundF3), cex.names = 0.5, las=3, ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRoundF3)), col="blue", lwd=2)

#read expersion: avg number counts per gene
rowSums(countsTableRoundF1)
mean(rowSums(countsTableRoundF1)) 
# 11930.81
median(rowSums(countsTableRoundF1)) 
# 2226

rowSums(countsTableRoundF3)
mean(rowSums(countsTableRoundF3)) 
# 11220.54
median(rowSums(countsTableRoundF3)) 
# 2144

apply(countsTableRoundF1,2,mean) # 1= rows, 2= columns
apply(countsTableRoundF1,1,mean)
hist(apply(countsTableRoundF1,1,mean), xlim=c(0,1000), breaks=10000) #change lim and breaks to see histogram better
apply(countsTableRoundF3,2,mean) # 1= rows, 2= columns
apply(countsTableRoundF3,1,mean)
hist(apply(countsTableRoundF3,1,mean), xlim=c(0,1000), breaks=10000) #change lim and breaks to see histogram better

###DESeq stuff:
#step1, create DESeq object and define experimanetal design with tilda(~)
dds1 <- DESeqDataSetFromMatrix(countData = countsTableRoundF1, colData = condsF1, design = ~ line + environment + line:environment)
dim(dds1) # has same number genes = 24362, and samples =16

dds3 <- DESeqDataSetFromMatrix(countData = countsTableRoundF3, colData = condsF3, design = ~ line + environment + line:environment)
dim(dds3) # has same number genes = 25379, and samples =16

#Filter genes with too few reads, keep avg greater then 10 reads per sample (10*16=160)
dds1 <-dds1[rowSums(counts(dds1)) > 160]
dim(dds1) #[1] 24362    16
dds3 <-dds3[rowSums(counts(dds3)) > 160]
dim(dds3) #[1] 25396   16


# Run DESeq model test for DEG
dds1 <- DESeq(dds1)
dds3 <- DESeq(dds3)

#list results generated
resultsNames(dds1)
#[1] "Intercept"                  "line_combined_vs_ambient"   "environment_HH_vs_AA"      
#[4] "linecombined.environmentHH"
resultsNames(dds3)
#[1] "Intercept"                  "line_combined_vs_ambient"   "environment_HH_vs_AA"      
#[4] "linecombined.environmentHH"

#PCA to visualize global gene expression patterns
vsd1 <- vst(dds1, blind=FALSE)
vsd3 <- vst(dds3, blind=FALSE)

data1 <- plotPCA(vsd1, intgroup=c("line","environment"), returnData=TRUE)
percentVar1 <- round(100 * attr(data1,"percentVar"))
data3 <- plotPCA(vsd3, intgroup=c("line","environment"), returnData=TRUE)
percentVar3 <- round(100 * attr(data3,"percentVar"))

ggplot(data1, aes(PC1,PC2, color=environment, shape=line)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) +
  theme_minimal()
ggplot(data3, aes(PC1,PC2, color=environment, shape=line)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar3[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar3[2],"% variance")) +
  theme_minimal()
# What patterns do we see? Clustering by groups, line and environment. 
# What gene expression results do we expect for each factor, main effects and/or interactions?
resInteraction1 <- results(dds1, alpha=0.05)
resInteraction1 <- resInteraction1[order(resInteraction1$padj),]
head(resInteraction1)
# baseMean   log2FoldChange             lfcSE             stat
# <numeric>        <numeric>         <numeric>        <numeric>
#TRINITY_DN115950_c0_g1 2245.96833710016 4.05236657027432 0.358490093639215 11.3039848023042
summary(resInteraction1)
#out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2839, 12%
# LFC < 0 (down)     : 1053, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 473, 1.9%

resInteraction3 <- results(dds3, alpha=0.05)
resInteraction3 <- resInteraction3[order(resInteraction3$padj),]
head(resInteraction3)
# baseMean    log2FoldChange             lfcSE
# <numeric>         <numeric>         <numeric>
#   TRINITY_DN142181_c0_g4  384.448386219129  3.11803328353206 0.398309170262656
summary(resInteraction3)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 271, 1.1%
# LFC < 0 (down)     : 60, 0.24%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2451, 9.7%
###############################################################

###############################################################
#TEST FOR EFFECT OF ENVIRONMENT
dds <- DESeqDataSetFromMatrix(countData = countsTableRoundF3, colData = condsF3, 
                              design = ~ line + environment)
#Likelyhood ratio test = LRT, taking out effect in reduced so removed line so left with Envior
dds <- DESeq(dds, test="LRT", reduced=~line)
# List the results you've generated
resultsNames(dds)
#[1] "Intercept"                "line_combined_vs_ambient" "environment_HH_vs_AA"

# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05)
head(resEnv) #not orded yet, next line gives ordered padj
resEnv <- resEnv[order(resEnv$padj),]
head(resEnv)
# baseMean    log2FoldChange             lfcSE
# <numeric>         <numeric>         <numeric>
#   TRINITY_DN121599_c1_g1   728.154992622927 -6.02233134180935 0.472156374478784
summary(resEnv)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 513, 2%
# LFC < 0 (down)     : 315, 1.2%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 491, 1.9%
#removing NA from resEnv
resEnv <- resEnv[!is.na(resEnv$padj),]
dim(resEnv)
#[1] 24772     6
#pulling out row names of 
degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 
length(degsEnv)
#828

###############################################################
#TEST FOR EFFECT OF LINE
dds <- DESeqDataSetFromMatrix(countData = countsTableRoundF3, colData = condsF3, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)
# [1] "Intercept"                "environment_HH_vs_AA"     "line_combined_vs_ambient"

resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)
# baseMean    log2FoldChange             lfcSE             stat
# <numeric>         <numeric>         <numeric>        <numeric>
#   TRINITY_DN132194_c0_g1 1229.94908215473   1.5390060139067 0.154824791799352  94.681394464729
summary(resLine)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 808, 3.2%
# LFC < 0 (down)     : 837, 3.3%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 981, 3.9
resLine <- resLine[!is.na(resLine$padj),]

degsline <- row.names(resLine[resLine$padj < 0.05,])

###############################################################
#TEST FOR INTERACTION
dds <- DESeqDataSetFromMatrix(countData = countsTableRoundF3, colData = condsF3, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)
# [1] "Intercept"                  "environment_HH_vs_AA"       "line_combined_vs_ambient"  
# [4] "environmentHH.linecombined"

resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)
# baseMean    log2FoldChange             lfcSE
# <numeric>         <numeric>         <numeric>
#   TRINITY_DN142181_c0_g4  384.448386219129  3.11803328353205 0.398309170262532
summary(resInt)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 235, 0.93%
# LFC < 0 (down)     : 48, 0.19%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2941, 12%
resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])

#Top gene from Env
d <-plotCounts(dds, gene="TRINITY_DN121599_c1_g1", intgroup = (c("line","environment")), returnData=TRUE)
d #prints data of gene
p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
#Top gene from Line
d <-plotCounts(dds, gene="TRINITY_DN132194_c0_g1", intgroup = (c("line","environment")), returnData=TRUE)
d #prints data of gene
p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
#Top gene from Interaction
d <-plotCounts(dds, gene="TRINITY_DN142181_c0_g4", intgroup = (c("line","environment")), returnData=TRUE)
d #prints data of gene
p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

# # By Environment
# library(pheatmap)
# topgenes <- head(rownames(resEnv),20)
# mat <- assay(vsd3)[topgenes,]
# mat <- mat - rowMeans(mat)
# df <- as.data.frame(colData(dds)[,c("line","environment")])
# pheatmap(mat, annotation_col=df)


library(eulerr)

# Total
length(degsEnv)  # 828
length(degsline)  # 1645
length(degsInt)  # 283

# Intersections
length(intersect(degsEnv,degsline))  # 141
length(intersect(degsEnv,degsInt))  # 14
length(intersect(degsInt,degsline))  # 32

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
828-14-141-7 # 666
1645-151-32-7 # 1455
283-14-32-7 # 230

#above did calculation for venn diagram: F3
fit1 <- euler(c("Env" = 666, "Line" = 1455, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))
fit1 <- euler(c("F1" = , "F2" = ))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

# Number unique
448-44-37-7 # 360
226-37-34-7 # 148
3854-44-34-7 # 3769

#above did calculation for venn diagram: F1
fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))
