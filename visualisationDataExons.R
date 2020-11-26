#### R script to visualise data from the workflow ####
  # Need the results from the R script "analyseData.R" 

# Remove all objects from workspace
rm(list=ls())


### Install Biocmanager and the depending packages ###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("pheatmap")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")

library("DESeq2")
library("pheatmap")
library("EnhancedVolcano")


### Exons data from the analysis ###

data_exons <- read.table("exons_analysis.txt", header = T, sep = "\t")
rownames(data_exons) = data_exons$row # Change the name of the rows


### Volcano Plot ###

EnhancedVolcano(data_exons,
                lab = data_exons$row,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,)

# Select exons significatively under or over transcript
listeExons <- data_exons$row[data_exons$pvalue<0.05]


### HeatMap ###

# Transform the dataframe into a matrix with the selection of the columns for the samples
count_matrix_exons <- as.matrix(data_exons[,c(8:15)]) 

# Metadata : needed to create DESeq object used for the pheatmap function
# Can be found at https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa
meta_Data <- read.table("SraRunTable.txt", sep="\t", header = T)
colData <- meta_Data[c("Run","LibraryLayout","sf3b1_mutation_status")]
# Create a new column indicating if each sample (in row) is a mutant ( SF3B1_mutated) or a wild-type (SF3B1_WT) for the SF3B1 gene.
colData$sf3b1_mutation_status_clean <- sapply(colData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),'SF3B1_WT','SF3B1_mutated')} )
colData$condition <- as.factor(colData$sf3b1_mutation_status_clean) # Type of information used in the DESeq funtion as a condition

# Combination of the exons data and the metadata into a DESeq object to be used in the DESeq function
dds <- DESeqDataSetFromMatrix(countData = count_matrix_exons,
                              colData = colData,
                              design = ~condition,
                              tidy = FALSE)

# Calculation with DESeq function --> results from the analyse of the data
dds <- DESeq(dds)
# Normalisation of the data (by default log2)
ntd <- normTransform(dds)
# Dataframe with information for the heatmap, for each sample, its status (mutated or wild type) and if the exons are paired. 
df <- as.data.frame(colData(dds)[,c("sf3b1_mutation_status_clean","LibraryLayout")])
# Heatmap 
pheatmap(assay(ntd)[listeExons,], cluster_rows=F, show_rownames=F,
         cluster_cols=F, annotation_col=df)


### Search of the genes' exons of the articles in our data from the workflow ###

# One gene of the article not found on the Internet (UQCC)

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
df_F8$gene <- c("F8")

exonsGAS8 <- exonsEtude$GAS8[is.na(exonsEtude$GAS8) == F]
df_GAS8 <- data_exons[data_exons$row %in% exonsGAS8,]
df_GAS8$gene <- c("GAS8")


# Combination of the genes' exons of the article with the exons available in our data
dfExonsArticleDonnees <- rbind(df_ABCC5, df_ADAM12,df_GUSBP11, df_ANKHD1, df_CRNDE, df_F8, df_GAS8)
write.table(dfExonsArticleDonnees, file = "data_exonsArticle.txt", quote = F, sep = "\t", row.names = F)

### HeatMap ### 
# Heatmap for the genes' exons found in the article

pheatmap(assay(ntd)[dfExonsArticleDonnees$row,], cluster_rows=F, show_rownames=T,
         cluster_cols=F, annotation_col=df)
