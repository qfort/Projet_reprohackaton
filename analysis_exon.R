

rm(list=ls()) ##Nettoyage de l'environnement

#### Installation des packages ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("pheatmap")
BiocManager::install("DESeq2")
library("DESeq2")
library("pheatmap")


#### ANALYSE POUR TOUS LES EXONS ####
#### Nettoyage des données ####

## Données sur les exons
countData <- read.table("exon_output_wt0.counts", sep ="\t", skip=1, header=T, row.names = 1)
colnames(countData) <- sub(".bam", "", colnames(countData))# renome les colonnes
colnames(countData) <- sub("bam_folder.", "", colnames(countData))#renome les colonnes
#et transforme en matrice.

## Metadonnées
meta_Data <- read.table("SraRunTable.txt", sep="\t", header = T)
colData <- meta_Data[c("Run","LibraryLayout","sf3b1_mutation_status")]
#creation d'une nouvelle colonne avec le status des tumeurs WT ou mutant de facon propre.
colData$sf3b1_mutation_status_clean <- sapply(colData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),'WT','mutated')} )

#### analyse des données
count_matrix <- as.matrix(countData[-(1:5)])#selectionne uniquement les colonnes avec les donnees de compte d'exon,

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ sf3b1_mutation_status_clean)#creation de l'objet utilisable pour l'analyse par deseq (et de l'utiliser par heatmap)

select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:400]#creation d'un objet permettant de selectionner les lignes a tracer dans la heatmap.
#On ordonne les colonne et affiches les 400 premiere.

df <- as.data.frame(colData(dds)[,c("sf3b1_mutation_status_clean","LibraryLayout")])#dataframe avec en correspondance l echantillon, son status SF3B1 et son type paired.
ntd <- normTransform(dds)#normalisation des donne (par default log2)

#### Affichage de la heatmap ####
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


#### ANALYSE POUR LE GENE ABCC5 ####
#selection des exons correspondant au gene ABCC5
data_table_ABCC5 <- subset(countData, countData$Chr=='3')#selection des exons du chromosome 3
#selection des exons selon leur position sur le chromosome 3
data_table_ABCC5<- subset(data_table_ABCC5,data_table_ABCC5$Start>='183984894')
data_table_ABCC5<- subset(data_table_ABCC5,data_table_ABCC5$End<='184017900')
count_matrix_ABCC5 <- as.matrix(data_table_ABCC5[-(1:5)])

dds_ABCC5 <- DESeqDataSetFromMatrix(countData = count_matrix_ABCC5,
                                    colData = colData,
                                    design = ~ sf3b1_mutation_status_clean)
select_ABCC5 <- order(rowMeans(counts(dds_ABCC5,normalized=FALSE)),
                      decreasing=TRUE)[1:6]
df_ABCC5 <- as.data.frame(colData(dds_ABCC5)[,c("sf3b1_mutation_status_clean","LibraryLayout")])
ntd_ABCC5 <- normTransform(dds_ABCC5) 
pheatmap(assay(ntd_ABCC5)[select_ABCC5,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df_ABCC5)