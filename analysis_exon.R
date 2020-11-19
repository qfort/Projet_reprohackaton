

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
dds <- DESeq(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]#creation d'un objet permettant de selectionner les lignes a tracer dans la heatmap.
#On ordonne les colonne et affiches les 400 premiere.

df <- as.data.frame(colData(dds)[,c("sf3b1_mutation_status_clean","LibraryLayout")])#dataframe avec en correspondance l echantillon, son status SF3B1 et son type paired.
ntd <- normTransform(dds)#normalisation des donne (par default log2)

#### Affichage de la heatmap ####
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)


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
dds_ABCC5 <- DESeq(dds_ABCC5)
select_ABCC5 <- order(rowMeans(counts(dds_ABCC5,normalized=FALSE)),
                      decreasing=TRUE)[1:6]
df_ABCC5 <- as.data.frame(colData(dds_ABCC5)[,c("sf3b1_mutation_status_clean","LibraryLayout")])
pheatmap(assay(ntd_ABCC5)[select_ABCC5,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df_ABCC5)


### ANALYSE volcano ###
BiocManager::install("EnhancedVolcano")
BiocManager::install("foreign")
library("EnhancedVolcano")
countData <- read.table("exon_output_wt0.counts", sep ="\t", skip=1, header=T, row.names = 1) 
countData <- countData[-(1:5)] ## éliminer les 5 premières colonnes inutiles

## Changer le nom des colonnes de CountData pour enlever le .bam
colnames(countData) <- sub(".bam", "", colnames(countData)) ## on enlève les 4 derniers caractères i.e le ".bam"
#rownames(countData)[duplicated(rownames(countData))]


## Récupération des métadonnées 
metaData <- read.csv('SraRunTable.txt', header = TRUE, sep = "\t")

metaData <- metaData[c("Run", "sf3b1_mutation_status")] #on ne garde que les colonnes d'intérêt

## CONSTITUTION DES DEUX GROUPES : SF3B1 muté et WT
sample_names <- metaData[,1] # récupération de l'index des échantillons 
SF3B1_WT <- c(grep("WT",metaData$sf3b1_mutation_status)) # on récupère l'index des échantillons WT pour SF3B1
SF3B1_mutated <- c(grep("WT",metaData$sf3b1_mutation_status, invert =T)) #on récupère le reste des index (avec le "invert" qui prend les index des échantillons ne répondant pas à la condition)

## construire une colonne condition dans metadata pour pouvoir établir le design
condition <- rep(c("0"),c(8))
condition[SF3B1_mutated] <- "SF3B1_mutated"
condition[SF3B1_WT] <- "SF3B1_WT"
condition <-  factor(condition,levels=c("SF3B1_mutated","SF3B1_WT"))

metaData$condition <- condition

## transformation de countData en matrice 
countData_mat <- data.matrix(countData)
## transformation des données en "DESeqDataSet objet"
dds <- DESeqDataSetFromMatrix(countData=countData_mat,
                              colData=metaData,
                              design=~ condition)


### Pre-filtering : enlever tous les gènes qui nous intéressent pas (faible nombre de counts)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

## DIFFERENTIAL EXPRESSION ANALYSIS
dds <- DESeq(dds) ## faire l'analyse diffé
res <- results(dds, contrast=c("condition","SF3B1_mutated","SF3B1_WT"), lfcThreshold = 1, alpha=0.05, tidy = T) ## obtenir le logFC en imposant l'ordre, tidy = T pour avoir les id dans une colonne
#res <- lfcShrink(dds, coef="SF3B1_mutated_vs_SF3B1_WT", type="apeglm") ## obtenir les log-FoldChange

## pour enlever les données NA 
res <- res[-which(is.na(res$log2FoldChange)),]
res <- res[-which(is.na(res$pvalue)),]



## VISUALISATION ###
## volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-4,
                FCcutoff = 2,)

## plotMA
#plotMA(dds) 

### Pour avoir un fichier excel en sortie 
## on prend les gènes les plus significatifs
#topgenes <- res[res[,"pvalue"] < 0.01, ]
# transformation en dataframe pour pouvoir filtrer
#topgenes <- as.data.frame(topgenes)

# on prend les gènes sur-régulés/exprimés
#topups <- topgenes %>% 
#  filter(log2FoldChange>0) %>% 
#  arrange(desc(log2FoldChange))

# on prend les gènes sous-régulés/exprimés
#topdowns <- topgenes %>% 
#  filter(log2FoldChange<0) %>% 
#   arrange(log2FoldChange)
