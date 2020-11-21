########### Etude de l'expression différentielle des gènes selon le groupe d'intérêt ##################


rm(list=ls()) ##Nettoyage de l'environnement

## Chargement des librairies à utiliser 
library("DESeq2")
library(tidyverse)
library(tibble)
library(EnhancedVolcano)
library(org.Hs.eg.db)

### To define the working directory
setwd('C:/Users/maryo/OneDrive/3A/COURS/AMI2B/Reprohackaton/Projet') 


#### PART1 : IMPORTING DATA ####
## Récupération de la table.counts
countData <- read.table("gene_output.counts", sep ="\t", skip=1, header=T, row.names = 1)
countData <- countData[-(1:5)] ## éliminer les 5 premières colonnes inutiles

## Changer le nom des colonnes de CountData pour enlever le .bam
colnames(countData) <- sub(".bam", "", colnames(countData)) ## on enlève les 4 derniers caractères i.e le ".bam"
#rownames(countData)[duplicated(rownames(countData))]


## To get metadata and keep the columns of interest
metaData <- read.csv('SraRunTable.txt', header = TRUE, sep = "\t")
metaData <- metaData[c("Run", "sf3b1_mutation_status")] #on ne garde que les colonnes d'intérêt


## CONSTITUTION DES DEUX GROUPES : SF3B1 muté et WT
## construire une colonne condition dans metadata pour pouvoir établir le design : chercher les échantillons qui sont WT
condition <- sapply(metaData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),"SF3B1_WT","SF3B1_mutated")}) #si pas WT -> SF3B1 mutated
# Transformer en factor pour pouvoir utiliser le design
metaData$condition <- factor(condition, levels = c("SF3B1_WT","SF3B1_mutated")) 


## Transformation de countData en matrice 
countData_mat <- data.matrix(countData)
## transformation des données en "DESeqDataSet objet"
dds <- DESeqDataSetFromMatrix(countData=countData_mat,
                              colData=metaData,
                              design=~ condition)




#### PART 2 : DIFFERENTIAL EXPRESSION ANALYSIS  ####
dds <- DESeq(dds) ## faire l'analyse différentielle des gènes 

res <- results(dds, tidy =T) ## obtenir le logFC en imposant l'ordre, tidy = T pour avoir les id dans une colonne

## pour enlever les données NA dans colonne logFoldChange et pvalue
res <- res[-which(is.na(res$log2FoldChange)),]
res <- res[-which(is.na(res$pvalue)),]

## volcano plot avec les ENSEMBL id
EnhancedVolcano(res,
                lab = res[,1],
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0,10),
                pCutoff = 0.1,
                FCcutoff = log2(1.5),
                labSize = 4.0)

### to keep only topgenes from our DE analysis
topgenes <- res[res[,"pvalue"] < 0.05, ]

topgenes <- as.data.frame(topgenes)
topups <- topgenes %>% 
  filter(log2FoldChange>(log2(1.5))) %>% 
  arrange(desc(log2FoldChange))
  
topdowns <- topgenes %>% 
  filter(log2FoldChange<(-(log2(1.5)))) %>% 
  arrange(log2FoldChange)

## ajouter les deux types de gènes ensemble dans un même fichier de sortie
fichier_sortie <- rbind(topups,topdowns)

## créer une colonne avec le nom des gènes
hs <- org.Hs.eg.db
ts <- select(hs, keys = fichier_sortie[,1], columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
ts <- ts[!duplicated(ts$ENSEMBL),]

## ajouter cette colonne dans le fichier de sortie
fichier_sortie$Gene_symbol <- ts$SYMBOL
colnames(fichier_sortie)[1] <- "ENSEMBL_id"

### écrire le fichier de sortie dans votre working directory
write.table(fichier_sortie, file="./DE_analysis.txt", sep = "\t", row.names = FALSE)






### PART3 : COMPARING WITH DE GENES FOUND IN THE ARTICLES ####
## récupération données de l'article
DE_article <- read.table('DE_articles.txt', header = TRUE, sep = "\t")

## récupération des gènes en commun dans une liste 
gene_list <- paste(fichier_sortie[,1], collapse="|")
common <- c(grep(gene_list, DE_article[,1], value =T))


## récupérer les gènes en commun dans leur analyse et les ordonner alphabétiquement pour avoir le même ordre de gènes
data_article <- DE_article[which(DE_article$EnsEMBLID %in% common),]
data_article <- data_article[order(data_article$EnsEMBLID),]

## To add a new column with the log2FC values
data_article$Fold.Change <- as.numeric(data_article$Fold.Change)
data_article$log2FoldChange <- sapply(data_article$Fold.Change, function(x) {-log2(x)})


## récupérer les gènes en commun dans notre analyse perso et les ordonner alphabétiquement pour avoir le même ordre de gènes
our_analysis <- fichier_sortie[which(fichier_sortie[,1] %in% common),]
our_analysis <- our_analysis[order(our_analysis$row),]

## fusion des deux dataframes
data_merge <- cbind(data_article, our_analysis$log2FoldChange, our_analysis$pvalue)
colnames(data_merge) <- c("ENSEMBL_id", "FAST.DB.STABLE_id", "Gene_Symbol", "Regulation", "Article.|FC|", 
                          "Article.pvalue", "Aliases&Synonyms", "Article.log2FC", "Our.log2FC", "Our.pvalue")
  
  
### écrire le fichier de sortie dans votre working directory
write.table(data_merge, file="./DE_genes_comparison.txt", sep = "\t", row.names = FALSE)
