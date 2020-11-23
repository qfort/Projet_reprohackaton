########### Differential expression analysis of genes according to the group of interest ##################

rm(list=ls()) ## environment cleaning

## Loading of the librairies 
library("DESeq2")
library(EnhancedVolcano)

### To define the working directory
setwd('/mnt/mydatalocal/') ## A MODIFIER !!!! 


#### PART1 : IMPORTING DATA ####
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

## Creation of a volcano plot with ENSEMBL id
EnhancedVolcano(res,
                lab = res[,1],
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(0,10),
                pCutoff = 0.05,
                FCcutoff = log2(1.5),
                labSize = 4.0)

### Keep only topgenes from our DE analysis
topgenes <- res[res[,"pvalue"] < 0.05, ]

### To convert topgenes into dataframe
topgenes <- as.data.frame(topgenes)
## to keep only the top up-regulated genes
topups <- topgenes[topgenes[,"log2FoldChange"]>(log2(1.5)),]
topups <- topups[order(-topups$log2FoldChange),]

## to keep only the top down-regulated genes
topdowns <- topgenes[topgenes[,"log2FoldChange"]<(-log2(1.5)),]
topdowns <- topdowns[order(topdowns$log2FoldChange),]

## Gather the two types of genes in a same output file
fichier_sortie <- rbind(topups,topdowns)

## Add this column in the output file
colnames(fichier_sortie)[1] <- "ENSEMBL_id"

### Write the output file in our working directory
write.table(fichier_sortie, file="./DE_analysis.txt", sep = "\t", row.names = FALSE)






### PART3 : COMPARING WITH DE GENES FOUND IN THE ARTICLES ####
## Collect data from the article
DE_article <- read.table('DE_articles.txt', header = TRUE, sep = "\t")

## Collect common genes in a list 
gene_list <- paste(fichier_sortie[,1], collapse="|")

gene_list2 <- paste(DE_article[,1], collapse="|")
gene_list2 <- gsub(" // ", "|", gene_list2)

common_art <- c(grep(gene_list, DE_article[,1], value =T))
common_our <- c(grep(gene_list2, fichier_sortie[,1], value =T))

## Collect the analysis of common genes from the article and order it alphabetically
data_article <- DE_article[which(DE_article$EnsEMBLID %in% common_art),]
data_article <- data_article[order(data_article$EnsEMBLID),]

## Add a new column with the log2FC values
data_article$Fold.Change <- as.numeric(data_article$Fold.Change)
data_article$log2FoldChange <- sapply(data_article$Fold.Change, function(x) {-log2(x)})


## Collect the common genes in our analysis and order it alphabetically to have the same order
our_analysis <- fichier_sortie[which(fichier_sortie[,1] %in% common_our),]
our_analysis <- our_analysis[order(our_analysis$ENSEMBL_id),]

## Fusion of the two dataframes
data_merge <- cbind(data_article, our_analysis$log2FoldChange, our_analysis$pvalue)
colnames(data_merge) <- c("ENSEMBL_id", "FAST.DB.STABLE_id", "Gene_Symbol", "Regulation", "Article.|FC|", 
                          "Article.pvalue", "Aliases&Synonyms", "Article.log2FC", "Our.log2FC", "Our.pvalue")
  
  
### Write the output file in our working directory
write.table(data_merge, file="./DE_genes_comparison.txt", sep = "\t", row.names = FALSE)
