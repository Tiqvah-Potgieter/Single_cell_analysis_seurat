#Step 1: Setup the Seurat Object

# download libraries 
install.packages("Seurat")
install.packages("dplyr")
install.packages("patchwork")
install.packages("ggplot2")
 #install.packages("SingleR")
 #install.packages("celldex")

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyverse)

#Load the data (sample1)

sample1 <- readRDS("/home/tpotgieter/cai_data/star_output/star/mtx_conversions/sample1_matrix.rds")

#look at the the structure of the rds file. 

str(sample1) #You will see an hierarchy of objects and data frames, which allows you to identify the components you need to extract for creating a Seurat object

object1 <- CreateSeuratObject(counts = sample1$RNA, project = "SeuratProject", min.cells = 3, 
                              min.features = 200) #note the $RNA 

object1                            # An object of class Seurat 
                                   #36601 features across 9285 samples within 1 assay 
                                   #Active assay: RNA (36601 features, 0 variable features)

#Step 2: Standard pre-processing workflow 
#Step 2.2.1: QC and selecting cells for further analysis

object1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")

#Show QC metrics for the first 5 cells

head(object1@meta.data, 5) 

#visualize the QC metrics as a violin plot 

VlnPlot(object1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)

plot1 <- FeatureScatter(object1, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position = "none")

plot2 <- FeatureScatter(object1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none")

plot1 + plot2

#filter away cells that have unique feature counts (genes) over 5,000 or less than 200. 

#Also filter away cells that have > 15% mitochondrial counts 

object1 <- subset(object1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

#Visualize QC metrics again after filtering cells 

VlnPlot(object1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)

plot1 <- FeatureScatter(object1, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position = "none")
plot2 <- FeatureScatter(object1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none")

plot1 + plot2 

#Step 3: Normalizing the data 

object1 <- NormalizeData(object1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

#set seed and put two plots in one figure 

set.seed(123)
par(mfrow=c(1,2))

#original expression distribution 

raw_geneExp = as.vector(object1[['RNA']]@counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp)

#expression distribution after normalization 

logNorm_geneExp = as.vector(object1[['RNA']]@data) %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
hist(logNorm_geneExp)

#Step 4: Identification of highly variable features 

object1 <- FindVariableFeatures(object1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

#Identify the 10 most highly variable genes 

top10 <- head(VariableFeatures(object1),10)

#plot variable features with and without labels 

plot1 <- VariableFeaturePlot(object1) + 
  theme(legend.position = "top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position = "none")

plot1 + plot2 

#Step 5: Scaling the data 

all.genes <- rownames(object1)
object1 <- ScaleData(object1, features = all.genes, verbose = FALSE)

#Step 6: Perform linear dimensional reduction 

object1 <- RunPCA(object1, features = VariableFeatures(object = object1), verbose = FALSE)

VizDimLoadings(object1, dims = 1:2, reduction = 'pca')

DimPlot(object1, reduction = "pca")

DimHeatmap(object1, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(object1, dims = 1:9, cells = 500, balanced = TRUE)

#Step 7: Cluster the cells 

object1 <- FindNeighbors(object1, dims = 1:20, verbose = FALSE)
object1 <- FindClusters(object1, resolution = 0.5, verbose = FALSE)

#look at the cluster ID of the first five cells 

head(Idents(object1), 5)

#Step 8: Run non-linear dimensional reduction (UMAP/tSNE)

object1 <- RunUMAP(object1, dims = 1:20, verbose = FALSE)

#Get the UMAP plot of the single cell clustering results

DimPlot(object1, reduction = "umap")

#Visualize it using a tSNE plot 

object1 <- RunTSNE(object1, dims = 1:20, verbose = FALSE)
DimPlot(object1, reduction = "tsne")

#Label individual clusters 

DimPlot(object1, reduction = "umap", label = TRUE)

plot <- DimPlot(object = object1)
LabelClusters(plot = plot, id = 'ident')

#Save the object 

saveRDS(object1, file = "/home/tpotgieter/Single_cell_analysis_seurat/sample1_processed.rds")
