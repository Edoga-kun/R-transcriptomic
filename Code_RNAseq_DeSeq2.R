# ================================================= #
# == RNAseq code by Edgar Figueredo Carmona ======= #
# ================================================= #


# ======================================================== #
# == 1. Read the matrix counts Data and the factor file == #
# ======================================================== #

countTable= read.table("Counts_matrix.txt", header=T,row.names=1)
head (countTable)
colData= read.table ("Factors.txt", header=T)
head (colData)

# ================================================================ #
# == 2. Filtering out low expressed genes (less than 32 counts) == #
# ================================================================ #

countTable1 = countTable [,-1]
countTable1[, "max"]= apply(countTable1[, 1:ncol(countTable1)], 1, max)
countTable1=countTable1[countTable1[,ncol(countTable1)]>32,]
countTable1= countTable1 [,-ncol(countTable1)]

# ============================================== #
# == 3. differential expression DEseq2 library == #
# ============================================== #

library(DESeq2)

dds<-DESeqDataSetFromMatrix(countData= countTable1,colData= colData,design= ~ Condition)

#This command line is to specify the reference group
dds$Condition <- relevel(dds$Condition, ref = "Control")

#####Calling differential expression
dds<-DESeq(dds)
res<-results(dds)
summary (res)

##### Save all the results of gene expression
write.table(as.data.frame(res),file="Differentially_expressed_genes.txt", sep = "\t", quote = F)

#####save the normalized counts per sample
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data_normalized_counts.txt", sep="\t", quote=F)



##### DE genes with log2FC >|1| and p-value< 0.05 (more lax filter)

res1= as.data.frame(res) #En este de aca sobre este archivo lo puedo filtrar sin necesidad del excel

DE_genes_pvalue= subset.data.frame(res1, padj<0.05 & abs(log2FoldChange)>1)
write.csv ( DE_genes_pvalue, file="results_DEseq_pval_FC1.csv" )



####Filter Normalized counts to get just DE genes

list= row.names(DE_genes_pvalue)
Filter <- subset(normalized_counts, rownames(normalized_counts) %in% list)
Filter= as.matrix(Filter)
write.table  (Filter, file="normalizedcounts_filter.txt")



# =========================== #
# == 4. TPM normalization  == #
# =========================== #

# Compute TPM for a read count matrix

tpm <- function(dfr,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#### 

counts= countTable[,2:5]

length =  countTable[,1]

TPM= tpm (counts, length)

write.table  (TPM, file="TPM.txt")


# ======================================== #
# == 7.Convert Ensembl-ID to GeneSymbol == #
# ======================================== #
library(biomaRt)
datasets <- listDatasets(useMart("ensembl"))
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
ensembl_genes <- rownames(Filter)
geneSymbols <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "description"),
  values= ensembl_genes,
  mart= mart)
colnames(geneSymbols) <- c("ensembl-ID","geneSymbol","description")

row.names (geneSymbols) = geneSymbols[,1]

test= merge(x = Filter , y = geneSymbols, by = "row.names", all = TRUE)
write.table( test, file="genes_names_matrix.txt", sep = "\t" )

# ====================================== #

