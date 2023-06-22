#Package installation and loading

install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")

devtools::install_github("jmw86069/splicejam")

#load packages 
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(data.table)
library(splicejam)
library(EnhancedVolcano)
library(DESeq2)

#setting work directory
setwd("D:/git_repositories/tag_seq5/")

#Importing transcript count data using both sample and sample1

dir = 'D:/git_repositories/tag_seq5/'

sampleData = paste0(dir, "sample1.csv")

sampleData = fread(sampleData)

rownames(sampleData) = sampleData$Samples
sampleData$Treatment = as.factor(sampleData$Treatment)
sampleData$Timepoint = as.factor(sampleData$Timepoint)

sampleData_VW5 = subset(sampleData, Treatment == 'VW5')


files = file.path(paste0(dir, "salmon"), list.files(paste0(dir, "salmon")), "quant.sf")
names(files) = list.files(paste0(dir, "salmon"))

gtf = "D:/git_repositories/tag_seq5/genome/mike_m_javanica_gtf"

tx2gene = makeTx2geneFromGtf(GTF = gtf, geneAttrNames = c("gene_id", "gene_name", "gene_type"), 
                             txAttrNames = c("transcript_id", "transcript_type"),
                             geneFeatureType = "gene",
                             txFeatureType = c("transcript", "mRNA"),
                             nrows = -1L,
                             verbose = FALSE)

#Importing salmon reads using tximport

txi = tximport(files, type = "salmon", tx2gene = tx2gene, txIn= TRUE, txOut= TRUE, countsFromAbundance = "no")


head(txi$counts)

#deseq2 analysis

dds = DESeqDataSetFromTximport(txi, colData = sampleData, ~Treatment)

head(dds)

dds = DESeq(dds)


cbind(resultsNames(dds))

rld = rlog(dds)
vsd = vst(dds)

pca_data <- plotPCA(vsd, intgroup = c("Treatment", "Timepoint"), returnData = TRUE)
percentVar <- round(100*attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color = group, shape = group)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC1: ", percentVar[1], "% variance")) + 
  coord_fixed()

res <- results(dds, name = "Treatment_VW5_0._vs_VW4_0.", alpha = 0.05)
summary(res)

write.csv(res, file='result_VW5_J2vsVW4_J2.csv')

down_VW5<- res[order(res$padj), ]
write.csv(down_VW5, file = "down_VW5_salmon_transcript.csv")

mcols(res)$description

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-9, 9),
                title = "VW5_J2 vs VW4_J2",
                subtitle = 'Differential Expression', 
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                pCutoff = 10e-3,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

resSig <- res[ which(res$padj<0.1), ]

head(resSig)
nrow(resSig)

write.csv(resSig, "significantDEgenes_VW4vsVW5.csv")


## Only VW5
dd_VW5 = DESeqDataSetFromTximport(txi, colData = sampleData_VW5, ~Timepoint)

