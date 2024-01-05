---
title: "SingleCellAnalysis_Alevin-Seurat"
format: html
editor: visual
---

# Analysis using Seurat

### Step 1: Setup the Seurat object

[Step 1.1: Downloading the packages]{.underline}

```{r}

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("patchwork")
#install.packages("Seurat")
```

[Step 1.2: Load the libraries]{.underline}

```{r}

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(tidyverse)
```

[Step 1.3: Reading in the data and initialize the Seurat object with the raw data]{.underline}

Note the format of the matrix files and load the data accordingly (different formats require different functions to load the the data). A .rds file uses the readRDS function.

```{r}

sample1 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample1_matrix.rds")
sample2 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample2_matrix.rds")
sample3 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample3_matrix.rds")
sample4 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample4_matrix.rds")
sample5 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample5_matrix.rds")
sample6 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample6_matrix.rds")
sample7 <- readRDS("/home/tpotgieter/cai_data/alevin_output/alevin/mtx_conversions/sample7_matrix.rds")


```

Look at the structure of the rds file

```{r}

str(sample1) 

#You will see an hierarchy of objects and data frames, which allows you to identify the components you need to extract for creating a Seurat object


```

[Step 1.4 : Create a Seurat object for each sample]{.underline}

[Step 1.4.1 : Create a function to create a Seurat object with varying samples as input]{.underline}

The CreateSeuratObject function in Seurat is used to create a Seurat object from a matrix of single-cell RNA-seq data.

The parameters for this object:

1.  count: This parameter should be a numeric matrix or a data frame where rows represent genes (features) and columns represent individual cells. Each element in the matrix should contain the expression value (e.g., read counts or UMI counts) for a specific gene in a specific cell.
2.  project: A character string that specifies the name of the project. This is used to label the Seurat object.
3.  min.cells: An integer specifying the minimum number of cells in which a gene must be detected to be retained. Genes not detected in at least this number of cells will be excluded from the Seurat object.
4.  min.features: An integer specifying the minimum number of features (genes) a cell must express to be retained. Cells expressing fewer features will be excluded from the Seurat object.

Note that there are more parameters but these ones are basic/standard ones

```{r}

#Create a function to create a Seurat object with varying samples as input
CreateSO <- function(sample_data, project_name) {
  
  # Extract the counts from the sample_data using $RNA
  count_input <- sample_data$RNA
  
  # Create the Seurat object 
  Seurat_object <- CreateSeuratObject(counts = count_input, project = project_name, min.cells = 3, min.features = 200)
  
  # Return the Seurat object 
  return(Seurat_object)
  
}


```

[Step 1.4.2: Create new Seurat objects using the function you created for each sample using a for loop]{.underline}

Note to use the assign function

The assign function takes two arguments:

-   The first argument is the name of the variable you want to create or modify (as a character string). In this case, it's the dynamically generated variable name created using paste.

-   The second argument is the value you want to assign to the variable.

```{r}

#Create new Seurat objects using the function you created for each sample using a for loop 

sample_data_list <- list(sample1, sample2, sample3, sample4, sample5, sample6, sample7)

#Iterate over the sample data and create Seurat object 
for (i in 1:length(sample_data_list)) { 
  
  project_name <- paste("sample", i, "project", sep = "")
  assign(paste0("object", i), CreateSO(sample_data_list[[i]], project_name))
  
  }
```

[Step 1.3.3: Look at the structure of the Seurat object]{.underline}

Features: There are 36,601 features. In the context of single-cell RNA-seq, "features" typically refer to genes or transcripts. These are the elements whose expression levels you've measured in your cells

Samples: There are 9,285 samples. In single-cell RNA-seq, these "samples" are individual cells. Each cell is treated as a sample, and the RNA expression data is captured for each of them

```{r}

object1
```

### Steps 2-4: Standard pre-processing workflow

The standard pre-processing workflow represents the selection and filtration of cells based on QC metrics, data normalization and scaling. A few QC metrics commonly used :

-   Low quality cells or empty droplets will often have very few genes

-   Cell doublets or multiplets may exhibit an aberrantly high gene count

-   The total number of molecules detected within a cell

-   Low quality/dying cells often exhibit extensive mitochondrial contamination

[Step 2.1: Set all genes starting with "MT-" as a set of mitochondrial genes]{.underline}

Calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features.

```{r}

object_data_list <- list(object1, object2, object3, object4, object5, object6, object7)

object_names <- c("object1", "object2", "object3", "object4", "object5", "object6", "object7")


for (i in 1:length(object_data_list)) { 
  
  object_new <- object_data_list[[i]]
  
  mt_percent <- PercentageFeatureSet(object_new, pattern = "^MT-")
  
 object_new[["percent.mt"]] <- mt_percent
 
 jab <-object_names[[i]]
 
 assign(jab,object_new)
 
  object_data_list[[i]] <- object_new
}


```

The number of unique genes and total molecules are automatically calculated during the "CreateSeuratObject()". They are stored in the object meta data.

```{r}

#show QC metrics 

head(object1@meta.data, 5)

```

[Step 2.2 : Visualize the QC metrics]{.underline}

nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell. And each dot in the following plots represents a cell.

```{r}

for (i in 1:length(object_data_list)) { 
  
vln_plot <- VlnPlot(object_data_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)

assign(paste0("vlnPlot", i), vln_plot)

}

vlnPlot1
```

We can use dot plots to show the relationship between nCount, nFeature and percent.mt

```{r}

 for (i in 1:length(object_data_list)) { 
   
   plot1 <- FeatureScatter(object_data_list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = "none")
   
   plot2 <- FeatureScatter(object_data_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = "none")
   
   assign(paste0("dotPlot", i), plot1 + plot2)
     
 }

dotPlot1
```

[Step 2.3: Apply filters to remove unwanted cells and visualize the QC metrics again]{.underline}

Filter away cells that have unique feature counts(genes) over 5 000 or less than 200. We also filter away cells that have \> 15% mitochondrial counts. This filtering criteria aims to remove low quality cells ( \>15% mitochondrial counts), potential doublets ( feature counts \> 5000) and potential empty droplets ( feature counts \< 200)

```{r}

#First filter cells 

for (i in 1:length(object_data_list)) { 
  
  object_data_list[[i]] <- subset(object_data_list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15) 
  
  assign(paste0("object", i), object_data_list[[i]])
  }


```

```{r}

#Visualize the QC metrics again after filtering cells 

for (i in 1:length(object_data_list)) { 
  
vln_plot <- VlnPlot(object_data_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)

assign(paste0("vlnPlot", i), vln_plot)

}

vlnPlot1

```

We can use dot plots to show the relationship between nCount, nFeature and percent.mt

```{r}

 for (i in 1:length(object_data_list)) { 
   
   plot1 <- FeatureScatter(object_data_list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = "none")
   
   plot2 <- FeatureScatter(object_data_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = "none")
   
   assign(paste0("dotPlot", i), plot1 + plot2)
     
 }

dotPlot1
```

[Step 3: Normalization]{.underline}

Normalization is a crucial step that aims to account for technical variability and differences in sequencing depth between individual cells, making the gene expression data more comparable and suitable for downstream analysis.

```{r}

 for (i in 1:length(object_data_list)) { 
   
 object_data_list[[i]] <- NormalizeData(object_data_list[[i]])
   
   assign(paste0("object", i), object_data_list[[i]])
     
 }
```

Normalization values are stored in object1\[\['RNA'\]\]\@data

```{r}

#set seed and put two plots in one figure 
set.seed(123)
par(mfrow = c(1,2))

# original expression distribution 
raw_geneExp = as.vector(object1[['RNA']]@counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp)

# expression distribution after normalization 

 logNorm_geneExp = as.vector(object1[['RNA']]@data) %>% sample(10000)
 logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
 hist(logNorm_geneExp)
```

[Step 4.1 Identification of highly variable features]{.underline}

We next calculate a subset of features that exhibit high cell-to-cell variation (they are highly expressed in some cells and lowly expressed in others). Focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

The FindVariableFeatures() function is used to identify genes that display high cell-to-cell variation. The 'vst' (Variance Stabilizing Transformation) parameter specifies the method used for feature selection. The nfeatures = 2000 determines the number of features/genes to retain, it specifies that the analysis should retain the top 2000 most variable genes.

```{r}

for (i in 1:length(object_data_list)) { 
  
  object_data_list[[i]] <- FindVariableFeatures(object_data_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  assign(paste0('object', i), object_data_list[[i]])
  
  # Identify the 10 most highly variable genes 
  top10 <- head(VariableFeatures(object_data_list[[i]]), 10) 
  
  # Plot variable features 
  
  plot1 <- VariableFeaturePlot(object_data_list[[i]]) + theme(legend.position = "top")
  
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.position = "none")
  
  assign(paste0('vfPlot', i), plot1 + plot2)
  
}

 vfPlot1
```

[Step 5.1: Scaling the data]{.underline}

This step shifts the expression of each gene, so that the mean expression across cells is 0. And scales the expression in each gene, so that the variance across cells is 1 ( this step gives equal weight in downstream analyses, so that highly expressed genes do not dominate). The results of this are stored in object\[\['RNA'\]\]\@scale.data

To make sure we don't leave any genes out of the heatmap later, we are scaling all genes and not only the highly variable genes

```{r}

for (i in 1:length(object_data_list)) { 
  
  all.genes <- rownames(object_data_list[[i]])
  
  object_data_list[[i]] <- ScaleData(object_data_list[[i]], features = all.genes, verbose = FALSE)
  
  
  }
```

[Step 6.1: Perform a linear dimensional reduction]{.underline}

Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using the features argument if you wish to choose a different subset.

Perform PCA on the scaled data.

```{r}

 for (i in 1:length(object_data_list)) { 
   
   object_data_list[[i]] <- RunPCA(object_data_list[[i]], features = VariableFeatures(object = object_data_list[[i]]), verbose = FALSE)
   
 }

 VizDimLoadings(object_data_list[[1]], dims = 1:2, reduction = "pca")
 DimPlot(object_data_list[[1]], reduction = "pca")
```

In particular DimHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses.

```{r}


DimHeatmap(object_data_list[[1]], dims = 1, cells = 500, balanced = TRUE) 
  
DimHeatmap(object_data_list[[1]], dims = 1:9, cells = 500, balanced = TRUE)

 
```

[Step 7.1: Cluster the cells]{.underline}

```{r}

for (i in 1:length(object_data_list)) { 
  
  object_data_list[[i]] <- FindNeighbors(object_data_list[[i]], dims = 1:20, verbose = FALSE)
  
  object_data_list[[i]] <- FindClusters(object_data_list[[i]], resolution = 0.5, verbose = FALSE)
  
  
  
  }
```

[Step 7.2: Run non-linear dimensional reduction (UMAP)]{.underline}

The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots.

```{r}

 for (i in 1:length(object_data_list)) { 
   
   object_data_list[[i]] <- RunUMAP(object_data_list[[i]], dims = 1:20, verbose = FALSE)
 }

 DimPlot(object_data_list[[1]], reduction = "umap")
 
 DimPlot(object_data_list[[1]], reduction = "umap", label = TRUE)
```
