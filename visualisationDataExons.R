##### Necessite les résultats du script analyseData.R ou du workflow associé #####


rm(list=ls()) ##Nettoyage de l'environnement

#### Installation des packages ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("pheatmap")
BiocManager::install("DESeq2")
library("DESeq2")
library("pheatmap")

BiocManager::install("EnhancedVolcano")
library("EnhancedVolcano")

#### Data 
data_exons <- read.table("exons_analysis.txt", header = T, sep = "\t")
rownames(data_exons) = data_exons$row

#### Volcano Plot 
EnhancedVolcano(data_exons,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 2,)

#Selection des exons significativements sur ou sous transcrits d'après le volcano plot
listeExons <- data_exons$row[data_exons$pvalue<0.05]


#### HeatMap
count_matrix_exons <- as.matrix(data_exons[,c(8:15)])
# Metadonnées
meta_Data <- read.table("SraRunTable.txt", sep="\t", header = T)
colData <- meta_Data[c("Run","LibraryLayout","sf3b1_mutation_status")]
# creation d'une nouvelle colonne avec le status des tumeurs WT ou mutant de facon propre.
colData$sf3b1_mutation_status_clean <- sapply(colData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),'SF3B1_WT','SF3B1_mutated')} )
colData$condition <- as.factor(colData$sf3b1_mutation_status_clean) #Pour le volcano plot

dds <- DESeqDataSetFromMatrix(countData = count_matrix_exons,
                              colData = colData,
                              design = ~condition,
                              tidy = FALSE)

#creation de l'objet utilisable pour l'analyse par deseq (et de l'utiliser par heatmap)
dds <- DESeq(dds)
ntd <- normTransform(dds)#normalisation des donne (par default log2)

df <- as.data.frame(colData(dds)[,c("sf3b1_mutation_status_clean","LibraryLayout")])#dataframe avec en correspondance l echantillon, son status SF3B1 et son type paired.

pheatmap(assay(ntd)[listeExons,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=TRUE, annotation_col=df)

#### Comparaison avec les exons des gènes relevés dans l'étude
# un gène de l'article QCC non retrouvé en cherchant sur internet

exonsEtude <- read.csv("exon_table_parsing.csv", header = T)
exonsEtude <- exonsEtude[,c(2:8)]
exonsEtude[exonsEtude == ""] <- NA

exonsABCC5 <- exonsEtude$ABCC5[is.na(exonsEtude$ABCC5) == F]
df_ABCC5 <- data_exons[data_exons$row %in% exonsABCC5,]
df_ABCC5$gene <- c("ABCC5")

exonsCRNDE <- exonsEtude$CRNDE[is.na(exonsEtude$CRNDE) == F]
df_CRNDE <- data_exons[data_exons$row %in% exonsCRNDE,]
df_CRNDE$gene <- c("CRNDE")

exonsGUSBP11 <- exonsEtude$GUSBP11[is.na(exonsEtude$GUSBP11) == F]
df_GUSBP11 <- data_exons[data_exons$row %in% exonsGUSBP11,]
df_GUSBP11$gene <- c("GUSBP11")

exonsANKHD1 <- exonsEtude$ANKHD1[is.na(exonsEtude$ANKHD1) == F]
df_ANKHD1 <- data_exons[data_exons$row %in% exonsANKHD1,]
df_ANKHD1$gene <- c("ANKHD1")

exonsADAM12 <- exonsEtude$ADAM12[is.na(exonsEtude$ADAM12) == F]
df_ADAM12 <- data_exons[data_exons$row %in% exonsADAM12,]
df_ADAM12$gene <- c("ADAM12")

exonsF8 <- exonsEtude$F8[is.na(exonsEtude$F8) == F]
df_F8 <- data_exons[data_exons$row %in% exonsF8,]
#df_F8$gene <- c("F8") ## --> 0 observations

exonsGAS8 <- exonsEtude$GAS8[is.na(exonsEtude$GAS8) == F]
df_GAS8 <- data_exons[data_exons$row %in% exonsGAS8,]
df_GAS8$gene <- c("GAS8")

#Combinaison des exons des gènes de l'article avec les exons de nos données issues du workflow
dfExonsArticleDonnees <- rbind(df_ABCC5, df_ADAM12, df_GUSBP11, df_ANKHD1, df_CRNDE, df_F8, df_GAS8)
write.table(dfExonsArticleDonnees, file = "data_exonsArticle.txt", quote = F, sep = "\t", row.names = F)

