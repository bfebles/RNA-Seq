---
title: "RNA-seq Analysis"
output: html_notebook
---

```{r}
if(!require('BiocManager',quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('DESeq2')
```

```{r}
library(DESeq2)
library(ggplot2)
library(tidyverse)
```

## Prepare data for DESeq2 object
```{r}
# import read counts data from featureCounts
countData<-read.csv('feature_count/featurecount.csv',header = TRUE, row.names = 1)

#view read count
head(countData)
dim(countData)
```
```{r}
#remove unwanted columns
countData <- countData%>%
  dplyr::select(-Chr,-Start,-End,-Strand,-Length)
```
```{r}
#filter out genes that has less than total of 50 mapped reads to avoid noise in further analysis
countData <- countData[rowSums(countData[,-1]) >=50,]
dim(countData)
#total of 18,483 genes
```

```{r}
# create dataframe with experiment labels with two conditions (IgG and Treatment)
condition=factor(c("IgG","IgG","IgG","IgG","IgG","IgG","IgG","IgG","IgG", "Treatment","Treatment","Treatment","Treatment","Treatment","Treatment","Treatment","Treatment"))

colData <- data.frame(row.names=colnames(countData),condition)
colData
```
## Create DESeq2 object
```{r}
#create input matrix
dds<-DESeqDataSetFromMatrix(countData , colData, ~condition)

#run DESeq
dds <- DESeq(dds)

```
## QC
```{r}
#plot fold change over the average expression level of all samples
plotMA(dds, ylim=c(-5,5))

#perform variance stabilizing transformation
vsdata <-vst(dds, blind = FALSE)

#plot PCA to check how samples cluster
plotPCA(vsdata, intgroup = 'condition', returnData = F)

# plot dispersion
plotDispEsts(dds)
```
## Differential gene expression
```{r}
#get differentially expressed genes
res <-results(dds, contrast = c('condition','IgG','Treatment'))
res

#order by adjusted p-value
resOrdered <-res[order(res$padj),]
resOrdered
```

```{r}
#remove NA and keep the padj value lower than 0.05 to keep only significant genes
sig<-resOrdered[!is.na(resOrdered$padj) &
resOrdered$padj<0.05 &
resOrdered$baseMean > 100 &
abs(resOrdered$log2FoldChange) > 1,]
sig

# create dataframe of sig
sig.df <- as.data.frame(sig)
write.csv(sig.df, file = 'feature_count/DESeq_results.csv')
```

```{r}
#load mouse annotation and ID library
BiocManager::install('org.Mm.eg.db')
library(org.Mm.eg.db)
```
```{r}
#add annotation column in the sig df
sig.df$Gene<-mapIds(org.Mm.eg.db, keys=rownames(sig.df), keytype = "ENSEMBL", column = "SYMBOL")
sig.df

#remove rows with NA
cbind(lapply(lapply(sig.df, is.na), sum))
sig.df<-na.omit(sig.df)
```
## Volcano plot 
```{r}
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
```
```{r}
vPlot<-EnhancedVolcano(sig.df, x = 'log2FoldChange', y='padj',lab = sig.df$Gene,colAlpha = 0.5,pointSize = 2,labSize = 4,
                title = 'Treatment vs IgG',
                legendPosition = 'right',legendLabSize = 8,legendIconSize = 5)
```
```{r}
png('vPlot.png', res = 300, width = 2000, height = 2000)
print(vPlot)
dev.off()
```

```{r}
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
```

```{r}
#get normalized counts from dds object
mat<-counts(dds, normalized = T)[rownames(sig.df),]
mat
#get z-score
mat.z <- t(apply(mat,1,scale))
colnames(mat.z) <-rownames(colData)
mat.z
```
## Heatmap
```{r}
# create dataframe of sig_h
Hmap<-Heatmap(mat.z, name = "z-score",cluster_rows = T, cluster_columns = T, row_labels = sig.df[rownames(mat.z),]$Gene, 
        column_labels = colnames(mat.z))
```
```{r}
png('heatmap.png', res = 300, width = 1500, height = 2500)
print(Hmap)
dev.off()
```
## Gene Ontology
```{r}
BiocManager::install('clusterProfiler')
BiocManager::install('AnnotationDbi')
```

```{r}
library(clusterProfiler)
library(AnnotationDbi)
```
```{r}
genes<-rownames(sig)
genes

GO_res <-enrichGO(gene=genes, OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL', ont = 'BP')
```
```{r}
GO.plot<-plot(barplot(GO_res, showCategory = 20, main = "Top 10 biological pathways", xlab="Gene Count"))
GO.plot
```
```{r}
png('GO.png', res = 300, width = 3000, height = 3500)
print(GO.plot)
dev.off()
```







