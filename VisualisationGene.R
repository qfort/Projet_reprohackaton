#### R script to visualise data from the workflow ####

# Remove all objects from workspace
rm(list=ls())


library("DESeq2")
library(EnhancedVolcano)
setwd('/mnt/mydatalocal') 

## Loading of table.counts
countData <- read.table("gene_output.counts", sep ="\t", skip=1, header=T, row.names = 1)
countData <- countData[-(1:5)] ## discard the 5 first columns to keep only the 8 samples in each column

## Change the name of the columns of CountData to remove .bam
colnames(countData) <- sub(".bam", "", colnames(countData)) ## discard the last 4 characters i.e ".bam"

## Get metadata and keep the columns of interest
metaData <- read.csv('SraRunTable.txt', header = TRUE, sep = "\t")
metaData <- metaData[c("Run", "sf3b1_mutation_status")] #we only keep columns of interest


## THE CONSTITUTION OF TWO GROUPS : SF3B1 mutated et WT
## Create a column for condition in metadata to establish the design : look for the mutated samples
condition <- sapply(metaData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),"SF3B1_WT","SF3B1_mutated")}) #if not WT -> SF3B1 mutated
# Turn into factor to use the design
metaData$condition <- factor(condition, levels = c("SF3B1_WT","SF3B1_mutated")) 


## Transformation of countData into matrix
countData_mat <- data.matrix(countData)
## Data transformation into "DESeqDataSet objet"
dds <- DESeqDataSetFromMatrix(countData=countData_mat,
                              colData=metaData,
                              design=~ condition)




#### PART 2 : DIFFERENTIAL EXPRESSION ANALYSIS  ####
dds <- DESeq(dds) ## to do the differential expression analysis of genes

res <- results(dds, tidy =T) ## tidy = T to have ENSEMBL id as the first column not as rownames 

## Remove the NA data of the logFoldChange and pvalue columns
res <- res[-which(is.na(res$log2FoldChange)),]
res <- res[-which(is.na(res$pvalue)),]


##### PART 3 : VOLCANO PLOT ##########
## volcano plot avec les ENSEMBL id
EnhancedVolcano(res,
                lab = res[,1],
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = log2(1.5),
                labSize = 4.0)


