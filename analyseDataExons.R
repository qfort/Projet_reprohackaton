rm(list=ls()) ##Nettoyage de l'environnement

#### Installation des packages ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("DESeq2")
library("DESeq2")

## Données sur les exons 
countData <- read.table("exon_output.counts", sep ="\t", skip=1, header=T, row.names = 1) 
countData$Exons <- rownames(countData) #Ajoute une colonne avec le nom des exons pour chaque ligne
colnames(countData) <- sub(".bam", "", colnames(countData))# renomme les colonnes
colnames(countData) <- sub("bam_folder.", "", colnames(countData))#renomme les colonnes
#et transforme en matrice.

## Metadonnées
meta_Data <- read.table("SraRunTable.txt", sep="\t", header = T)
colData <- meta_Data[c("Run","LibraryLayout","sf3b1_mutation_status")]
#creation d'une nouvelle colonne avec le status des tumeurs WT ou mutant de facon propre.
colData$sf3b1_mutation_status_clean <- sapply(colData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),'SF3B1_WT','SF3B1_mutated')} )
colData$condition <- as.factor(colData$sf3b1_mutation_status_clean) #Pour le volcano plot

####Analyse des données 
count_matrix <- as.matrix(countData)
count_matrix_exons <- as.matrix(countData[,c(6:13)])
dds <- DESeqDataSetFromMatrix(countData = count_matrix_exons,
                              colData = colData,
                              design = ~condition,
                              tidy = FALSE)
#creation de l'objet DESeq
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","SF3B1_mutated","SF3B1_WT"), lfcThreshold = 1, alpha=0.05, tidy = T) ## obtenir le logFC en imposant l'ordre, tidy = T pour avoir les id dans une colonne
#res <- lfcShrink(dds, coef="SF3B1_mutated_vs_SF3B1_WT", type="apeglm") ## obtenir les log-FoldChange
# pour enlever les données NA 
res <- res[-which(is.na(res$pvalue)),]


#Obtention du dataframe final avec la fusion des dataframes
data_output <- merge(res, countData[,c(6:14)], by.x = "row", by.y = "Exons", all = F)
write.table(data_output, file = "exons_analysis.txt", quote = F, sep = "\t", row.names = F)


