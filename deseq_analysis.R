## Setting up working directory

rm(list=ls())
setwd("~/git_repositories/Tag_seq_javanica/Input_files")


#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

BiocManager::install("vsn")

#install.packages("ggVennDiagram")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

BiocManager::install("DESeq2")


install.packages("BiocManager")

#BiocManager::install("rrvgo")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

## Getting all the libraries that we need
library(DESeq2)
library(RColorBrewer)
library(EnhancedVolcano)
library(vsn)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(gplots)
library(VennDiagram)
library(RColorBrewer)
#library(rrvgo)
library(topGO)
library(dplyr)
library(tibble)

# Getting the read count data

exp_data <- read.delim("only_javanica.txt", row.names = 1, sep = "\t")
#rownames(exp_data) <-gsub("transcript:", "",rownames(exp_data)) ## Removing the word transcript and replaces it with nothing
exp_data <- as.matrix(exp_data) #DESeq2 pipeline expects the count data as a matrix  
head(exp_data)
dim(exp_data)

# Constructing metadata with rows corresponding to the columns of the count data (in this case exp_data)

treatment <- c("VW4","VW4","VW4","VW4","VW4","VW4","VW4_J2","VW4_J2", "VW4_J2", "VW5","VW5","VW5", "VW5", "VW5", "VW5", "VW5_J2", "VW5_J2", "VW5_J2" )
timepoint <- c(28, 28, 28, 7, 7, 7, 0, 0, 0, 28, 28, 28, 7, 7, 7, 0, 0, 0 )

group <- factor(paste0(treatment, ".", timepoint))

metadata <- data.frame("group" = group, row.names= colnames(exp_data))
metadata

table(metadata$group)


# Making DEseq data format

dds <- DESeqDataSetFromMatrix(exp_data, metadata, ~group)
#calculate the linear correction factors for each sample:
dds <- estimateSizeFactors(dds)
sizeFactors(dds)## gives the correction factors for each sample 
normalized_dds <- counts(dds, normalized = TRUE)

# Data transformation and visualization

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
head(assay(vsd))

# PCA plot

pca_data <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE)
percentVar <- round(100*attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color = group, shape = group)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC1: ", percentVar[1], "% variance")) + 
  coord_fixed()

#Heat map

sampledistance <- dist(t(assay(vsd)))
SampledistanceMatrix <- as.matrix(sampledistance)
rownames(SampledistanceMatrix) <- paste(vsd$group)
colnames(SampledistanceMatrix) <- paste(vsd$group)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(SampledistanceMatrix, 
         clustering_distance_rows = sampledistance, 
         clustering_distance_cols = sampledistance, 
         col = colors)

## Differential expression analysis

dd1 <- DESeq(dds)
head(dd1)
dd1$group
########################################################################

## Comparing VW5_J2 against VW4_J2 (control)

res1 <- results(dd1, contrast = c("group","VW5_J2.0", "VW4_J2.0"), alpha = 0.05)
summary(res1)
#THERE ARE 134 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res1$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res1$pvalue))

EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW5_J2 vs VW4_J2',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)
#DESeq2 uses the so-called Benjamini-Hochberg (BH) adjustment; in brief, this method calculates for each gene an adjusted p value which answers the following question: if one called significant all genes with a p value less than or equal to this gene's p value threshold, what would be the fraction of false positives (the false discovery rate, FDR) among them (in the sense of the calculation outlined above)?
  
##These values, called the BH-adjusted p values, are given in the column padj of the results object.
##Hence, if we consider a fraction of 10% false positives acceptable, we can consider all genes with an adjusted p value below 10%=0.1 as significant. How many such genes are there?

sum(res1$padj<0.1, na.rm = TRUE)  #20

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig1 <- res1[ which(res1$padj<0.1), ]
write.csv(resSig1, file = "DEgenes_VW5J2vsVW4J2.csv")

nrow(resSig1) #30 significant genes

#only down regulated genes 

down_VW5_J2<- resSig1[order(resSig1$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW5_J2 <- subset(down_VW5_J2, log2FoldChange <=-1.5) #2 highly down regulated
head(down_VW5_J2)

#only up regulated genes 
up_VW5_J2 <- resSig1[order(resSig1$log2FoldChange, decreasing = TRUE), ]
up_VW5_J2 <- subset(up_VW5_J2, log2FoldChange >=1.5) #8 highly up regulated

##############################################################################################################################################

#Comparing VW5_7days with VW4_J2 

res2 <- results(dd1, contrast = c("group","VW5.7", "VW4_J2.0"), alpha = 0.05)
summary(res2)

#THERE ARE 575 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res2$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res2$pvalue))

EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW5_7days vs VW4_J2(control)',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

sum(res1$padj<0.1, na.rm = TRUE)  #20

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig2 <- res2[ which(res2$padj<0.1), ]
nrow(resSig2) #404 significant genes
write.csv(resSig2, file = "DEgenes_VW5_7vsVW4J2.csv")


#only down regulated genes 

down_VW5_7<- resSig2[order(resSig2$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW5_7 <- subset(down_VW5_7, log2FoldChange <=-1.5) #166 highly down regulated
head(down_VW5_7)

#only up regulated genes 
up_VW5_7 <- resSig2[order(resSig2$log2FoldChange, decreasing = TRUE), ]
up_VW5_7 <- subset(up_VW5_7, log2FoldChange >=1.5) #181 highly up regulated

###############################################################################################################################################

#Comparing VW5_28 days vs VW4_J2

res3 <- results(dd1, contrast = c("group","VW5.28", "VW4_J2.0"), alpha = 0.05)
summary(res3)

#THERE ARE 1148 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res3$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res3$pvalue))

EnhancedVolcano(res3,
                lab = rownames(res3),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW5_28days vs VW4_J2(control)',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

sum(res3$padj<0.1, na.rm = TRUE)  #1099

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig3 <- res3[ which(res3$padj<0.1), ]
nrow(resSig3) #404 significant genes
write.csv(resSig3, file = "DEgenes_VW5_28vsVW4J2.csv")


#only down regulated genes 

down_VW5_28<- resSig3[order(resSig3$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW5_28 <- subset(down_VW5_28, log2FoldChange <=-1.5) #525 highly down regulated
nrow(down_VW5_28)

#only up regulated genes 
up_VW5_28 <- resSig3[order(resSig3$log2FoldChange, decreasing = TRUE), ]
up_VW5_28 <- subset(up_VW5_28, log2FoldChange >=1.5) #164 highly up regulated
nrow(up_VW5_28)
##############################################################################################################################################################

#Comparing VW4_7days with VW4_J2 

res4 <- results(dd1, contrast = c("group","VW4.7", "VW4_J2.0"), alpha = 0.05)
summary(res4)

#THERE ARE 527 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res4$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res4$pvalue))

EnhancedVolcano(res4,
                lab = rownames(res4),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW4_7days vs VW4_J2(control)',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

sum(res4$padj<0.1, na.rm = TRUE)  #366

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig4 <- res4[ which(res4$padj<0.1), ]
nrow(resSig4) # significant genes
write.csv(resSig4, file = "DEgenes_VW4_7vsVW4J2.csv")


#only down regulated genes 

down_VW4_7<- resSig4[order(resSig4$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW4_7 <- subset(down_VW4_7, log2FoldChange <=-1.5) #145 highly down regulated
head(down_VW4_7)
nrow(down_VW4_7)

#only up regulated genes 
up_VW4_7 <- resSig4[order(resSig4$log2FoldChange, decreasing = TRUE), ]
up_VW4_7 <- subset(up_VW4_7, log2FoldChange >=1.5) #162 highly up regulated
nrow(up_VW4_7)

#############################################################################################################################################################

res5 <- results(dd1, contrast = c("group","VW4.28", "VW4_J2.0"), alpha = 0.05)
summary(res4)

#THERE ARE 1128 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res5$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res5$pvalue))

EnhancedVolcano(res5,
                lab = rownames(res5),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW4_28days vs VW4_J2(control)',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

sum(res5$padj<0.1, na.rm = TRUE)  #1079 significant genes

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig5 <- res5[ which(res5$padj<0.1), ]
nrow(resSig5) # significant genes
write.csv(resSig5, file = "DEgenes_VW4_28vsVW4J2.csv")


#only down regulated genes 

down_VW4_28<- resSig5[order(resSig5$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW4_28 <- subset(down_VW4_28, log2FoldChange <=-1.5) #518 highly down regulated
head(down_VW4_28)
nrow(down_VW4_28)

#only up regulated genes 
up_VW4_28 <- resSig5[order(resSig5$log2FoldChange, decreasing = TRUE), ]
up_VW4_28 <- subset(up_VW4_28, log2FoldChange >=1.5) #169 highly up regulated
nrow(up_VW4_28)

##############################################################################################################################################################

res6 <- results(dd1, contrast = c("group","VW5.28", "VW5_J2.0"), alpha = 0.05)
summary(res6)

#THERE ARE 1106 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res6$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res6$pvalue))

EnhancedVolcano(res6,
                lab = rownames(res6),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW5_28days vs VW5_J2(control)',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

sum(res6$padj<0.1, na.rm = TRUE)  #1053 significant genes

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig6 <- res6[ which(res6$padj<0.1), ]
nrow(resSig6) #1053 significant genes
write.csv(resSig6, file = "DEgenes_VW5_28vsVW5J2.csv")


#only down regulated genes 

down_VW5_28_VW5<- resSig6[order(resSig6$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW5_28_VW5 <- subset(down_VW5_28_VW5, log2FoldChange <=-1.5) #479 highly down regulated
head(down_VW5_28_VW5)
nrow(down_VW5_28_VW5)

#only up regulated genes 
up_VW5_28_VW5 <- resSig6[order(resSig6$log2FoldChange, decreasing = TRUE), ]
up_VW5_28_VW5 <- subset(up_VW5_28_VW5, log2FoldChange >=1.5) #178 highly up regulated
nrow(up_VW5_28_VW5)

##############################################################################################################################################################

res7 <- results(dd1, contrast = c("group","VW5.7", "VW5_J2.0"), alpha = 0.05)
summary(res7)

#THERE ARE 561 GENES WITH pvalue below 0.05 among 7347 genes for which the test succeeded in reporting a p-value:
sum(res7$pvalue < 0.05 , na.rm = TRUE)
table(is.na(res7$pvalue))

EnhancedVolcano(res7,
                lab = rownames(res7),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = 'VW5_7days vs VW5_J2(control)',
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-5,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

sum(res7$padj<0.1, na.rm = TRUE)  #416 significant genes

#Subseting the results table to these genes and then sort by log2 fold chage estimate to get the significant genes with the strongest down-regulation and up-regulation

resSig7 <- res7[ which(res7$padj<0.1), ]
nrow(resSig7) # 416 significant genes
write.csv(resSig7, file = "DEgenes_VW5_7daysvsVW5J2.csv")


#only down regulated genes 

down_VW5_7_VW5<- resSig7[order(resSig7$log2FoldChange), ]

#sub-setting by log2foldchange <=-1.5

down_VW5_7_VW5 <- subset(down_VW5_7_VW5, log2FoldChange <=-1.5) #164 highly down regulated
head(down_VW5_7_VW5)
nrow(down_VW5_7_VW5)

#only up regulated genes 
up_VW5_7_VW5<- resSig7[order(resSig7$log2FoldChange, decreasing = TRUE), ]
up_VW5_7_VW5 <- subset(up_VW5_7_VW5, log2FoldChange >=1.5) #185 highly up regulated
nrow(up_VW5_7_VW5)

##############################################################################################################################################################
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

## Venn Diagram to compare the differentially expressed genes 

#sig_VW5J2 <- read.csv('DEgenes_VW5J2vsVW4J2.csv')
#sig_VW5_28 <- read.csv("DEgenes_VW5_28vsVW4J2.csv")
#sig_VW5_7 <- read.csv("DEgenes_VW5_7vsVW4J2.csv")

sig_VW4_28 <- read.csv("DEgenes_VW4_28vsVW4J2.csv")
nrow(sig_VW4_28)
sig_VW4_7 <- read.csv('DEgenes_VW4_7vsVW4J2.csv')
nrow(sig_VW4_7)
sig_VW5_7_VW5 <- read.csv("DEgenes_VW5_7daysvsVW5J2.csv")
nrow(sig_VW5_7_VW5)
sig_VW5_28_VW5 <- read.csv("DEgenes_VW5_28vsVW5J2.csv")
nrow(sig_VW5_28_VW5)

#set1 <-sig_VW5J2$X
#set2 <-sig_VW5_28$X
#set3<- sig_VW5_7$X
set4<- sig_VW4_28$X
set5<- sig_VW4_7$X
set6<- sig_VW5_28_VW5$X
set7<- sig_VW5_7_VW5$X


display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n =8, name ="Set2" )

display_venn(x = list(set5, set7), category.names = c("7dpi_VW4", "7dpi_VW5"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("#FC8D62" ,"#8DA0CB"))

display_venn(x = list(set4, set6), category.names = c("28dpi_VW4", "28dpi_VW5"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("#66C2A5" ,"#E78AC3"))

display_venn(x = list(set5, set4), category.names = c("7dpi_VW5", "28dpi_VW5"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("#FC8D62" ,"#8DA0CB"))

display_venn(x = list(set7, set6), category.names = c("7dpi_VW4", "28dpi_VW4"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = c("#FC8D62" ,"#8DA0CB"))

