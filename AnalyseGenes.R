########### Study of the differential expression of genes according to the group of interest ##################


rm(list=ls()) ## cleaning 

## Loading of the librairies 
library("DESeq2")
library(tidyverse)
library(tibble)
library(EnhancedVolcano)
library(org.Hs.eg.db)

### To define the working directory
setwd('C:/Users/maryo/OneDrive/3A/COURS/AMI2B/Reprohackaton/Projet') 


#### PART1 : IMPORTING DATA ####
## Loading of table.counts
countData <- read.table("gene_output.counts", sep ="\t", skip=1, header=T, row.names = 1)
countData <- countData[-(1:5)] ## éliminer les 5 premières colonnes inutiles

## Change the name of the columns of CountData to remove .bam
colnames(countData) <- sub(".bam", "", colnames(countData)) ## on enlève les 4 derniers caractères i.e le ".bam"
#rownames(countData)[duplicated(rownames(countData))]


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

res <- results(dds, tidy =T) ## to getlogFC by imposing order, tidy = T to have the ID in the column

## Remove the NA data of the logFoldChange and pvalue columns
res <- res[-which(is.na(res$log2FoldChange)),]
res <- res[-which(is.na(res$pvalue)),]

## Creation of a volcano plot with ENSEMBL id
EnhancedVolcano(res,
                lab = res[,1],
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,10),
                pCutoff = 0.1,
                FCcutoff = log2(1.5),
                labSize = 4.0)

### Keep only topgenes from our DE analysis
topgenes <- res[res[,"pvalue"] < 0.05, ]

topgenes <- as.data.frame(topgenes)
topups <- topgenes %>% 
  filter(log2FoldChange>(log2(1.5))) %>% 
  arrange(desc(log2FoldChange))
  
topdowns <- topgenes %>% 
  filter(log2FoldChange<(-(log2(1.5)))) %>% 
  arrange(log2FoldChange)

## Gather the two types of genes in a same output file
fichier_sortie <- rbind(topups,topdowns)

## Create a column with the name of the genes
hs <- org.Hs.eg.db
ts <- select(hs, keys = fichier_sortie[,1], columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
ts <- ts[!duplicated(ts$ENSEMBL),]

## Add this column in the output file
fichier_sortie$Gene_symbol <- ts$SYMBOL
colnames(fichier_sortie)[1] <- "ENSEMBL_id"

### Write the output file in our working directory
write.table(fichier_sortie, file="./DE_analysis.txt", sep = "\t", row.names = FALSE)






### PART3 : COMPARING WITH DE GENES FOUND IN THE ARTICLES ####
## Collect data from the article
DE_article <- read.table('DE_articles.txt', header = TRUE, sep = "\t")

## Collect common genes in a list 
gene_list <- paste(fichier_sortie[,1], collapse="|")
common <- c(grep(gene_list, DE_article[,1], value =T))


## Collect the analysis of common genes from the article and order it alphabetically
data_article <- DE_article[which(DE_article$EnsEMBLID %in% common),]
data_article <- data_article[order(data_article$EnsEMBLID),]

## Add a new column with the log2FC values
data_article$Fold.Change <- as.numeric(data_article$Fold.Change)
data_article$log2FoldChange <- sapply(data_article$Fold.Change, function(x) {-log2(x)})


## Collect the common genes in our analysis and order it alphabetically to have the same order
our_analysis <- fichier_sortie[which(fichier_sortie[,1] %in% common),]
our_analysis <- our_analysis[order(our_analysis$row),]

## Fusion of the two dataframes
data_merge <- cbind(data_article, our_analysis$log2FoldChange, our_analysis$pvalue)
colnames(data_merge) <- c("ENSEMBL_id", "FAST.DB.STABLE_id", "Gene_Symbol", "Regulation", "Article.|FC|", 
                          "Article.pvalue", "Aliases&Synonyms", "Article.log2FC", "Our.log2FC", "Our.pvalue")
  
  
### Write the output file in our working directory
write.table(data_merge, file="./DE_genes_comparison.txt", sep = "\t", row.names = FALSE)
