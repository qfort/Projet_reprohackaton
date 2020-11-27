#### R script to analyze data from the workflow ####

# Remove all objects from workspace
rm(list=ls())


### Install Biocmanager and the depending packages ###

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#BiocManager::install("DESeq2")
library("DESeq2")

setwd('.')
### Exons data ###

# Import exons data from the workflow : output file named "exon_output.counts"
countData <- read.table("exon_output.counts", sep ="\t", skip=1, header=T, row.names = 1) 
countData$Exons <- rownames(countData) # Add a column with the exons' names for later analysis' purposes
colnames(countData) <- sub(".bam", "", colnames(countData)) # Rename columns
colnames(countData) <- sub("bam_folder.", "", colnames(countData))


### Metadata ###

# Import metadata found at https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa
# Information on each studied samples and the type of data 
meta_Data <- read.table("SraRunTable.txt", sep="\t", header = T)
meta_Data <- meta_Data[c(1,2,3,4,5,6,8,7),] # Rows for SRR628589 and SRR628588 are not in a good order 

colData <- meta_Data[c("Run","LibraryLayout","sf3b1_mutation_status")] # Extract the columns of interest for our analysis
# Create a new column indicating if each sample (in row) is a mutant ( SF3B1_mutated) or a wild-type (SF3B1_WT) for the SF3B1 gene.
colData$sf3b1_mutation_status_clean <- sapply(colData$sf3b1_mutation_status, function(x){ifelse(grepl(pattern = "WT",x),'SF3B1_WT','SF3B1_mutated')} )
colData$condition <- as.factor(colData$sf3b1_mutation_status_clean) # Transform in a new column the type of the information about the mutation status : it will be used for the DESeq function

### Data analysis ###
count_matrix <- as.matrix(countData) # Transform the dataframe into a matrix so it can be used in the DESeq functions
count_matrix_exons <- as.matrix(countData[,c(6:13)]) # Transform the dataframe into a matrix so it can be used in the DESeq functions. 
# Selection of the columns with the information for each sample (8 in total), the other columns are not useful for our analysis


# Combination of the exons data and the metadata into a DESeq object to be used in the DESeq function
dds <- DESeqDataSetFromMatrix(countData = count_matrix_exons,
                              colData = colData,
                              design = ~condition,
                              tidy = FALSE)

# Calculation with DESeq function --> results from the analyse of the data
dds <- DESeq(dds)

# Extraction of the information contained in the DESeq object in a usable dataframe 
res <- results(dds, contrast=c("condition","SF3B1_mutated","SF3B1_WT"), tidy = T) 

# Remove N/A from the pvalue 
res <- res[-which(is.na(res$pvalue)),]

# Final dataframe obtained from the merging of our two main dataframes : different information to be visualised in the visualisationDataExons.R script
data_output <- merge(res, countData[,c(6:14)], by.x = "row", by.y = "Exons", all = F)
write.table(data_output, file = "exons_analysis.txt", quote = F, sep = "\t", row.names = F)


