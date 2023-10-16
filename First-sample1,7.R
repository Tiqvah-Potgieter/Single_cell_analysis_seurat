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






