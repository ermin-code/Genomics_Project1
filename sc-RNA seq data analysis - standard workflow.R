#script to perform standard workflow steps to analyze single cell RNA-Seq data
#data: 10 Human PBMCs
#data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728

setwd('C:/Users/ermin/Documents/Bioinformatics Shift Project/Website Portfolio Project/Bioinformatics Portfolio Website/genomics/genomics_project1/GSE150728')

#loading libraries
library(Seurat)
library(tidyverse)

#PERFORMING STANDARD WORKFLOW ON sc-RNA DATA FROM PBMCS OF A HEALTHY PATIENT LABELED 'H1'

#loading the NSCLC dataset from downloaded RDS file
h1pbmcs.sparse.m <- readRDS('C:/Users/ermin/Documents/Bioinformatics Shift Project/Website Portfolio Project/Bioinformatics Portfolio Website/genomics/genomics_project1/GSE150728/Raw Data/GSE150728_RAW/H1/GSM4557334_HIP002_cell.counts.matrices.rds')

#running str command to see all fields in the class
str(h1pbmcs.sparse.m)

#selecting gene expression data and storing it into variable 'h1.counts.exon'
h1pbmcs.counts.exon <- h1pbmcs.sparse.m$exon

#display first 10 rows of the h1.counts.exon matrix
head(h1pbmcs.counts.exon, n=10)

#initializing Seurat object with the raw data
h1pbmcs.seurat.obj <- CreateSeuratObject(counts = h1pbmcs.counts.exon, project = "H-PBMCs", min.cells = 3, min.features = 200)

#running str command to see all fields in the class Note: Check out previous video to see what information is stored in Seurat object
str(h1pbmcs.seurat.obj)

#22575 features across 5150 samples

# 1. QUALITY CONTROL 
# Since this is a raw matrix and it hasn't been filtered, we want to filter out low quality cells

# Viewing Seurat object meta data

View(h1pbmcs.seurat.obj@meta.data)

  # A. % Mitochondrial Genes - low quality cells have higher mitochondrial gene contamination

    h1pbmcs.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(h1pbmcs.seurat.obj, pattern = "^MT-")

    # Viewing meta data that should include % of mitochondrial genes column
    
    View(h1pbmcs.seurat.obj@meta.data)
    
    # Visualizing the meta data using Violin plot. The plot includes three columns with number of features, number of molecules and percent of mitochondrial genes
    
    VlnPlot(h1pbmcs.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # Visualizing number of RNA molecules vs number of features via Scatter  - majority of good quality cells should follow straight line
    
    FeatureScatter(h1pbmcs.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
 
  # B. Filtering low quality cells based on the number of genes and cells that are having higher mitochondrial percentages Note: not included is checking doublets using doublet finder. Search videos on doublet finder.
    
      #B1. Filtering data based on cells that having greater than 200 genes and less than 2500 and keeping only those cells that have mitochondrial percentage less than 5%
    
      h1pbmcs.seurat.obj <- subset(h1pbmcs.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)
      
      #B2. Normalizing the data in order to compare gene expression across multiple cells
      
      # We divide gene expression measurements in each cell by the total expression and then multiply it by a scaling factor and then log transform it.
      
      #h1pbmcs.seurat.obj <-NormalizeData(h1pbmcs.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
      
      #Using normalization with default parameters
      
      h1pbmcs.seurat.obj <- NormalizeData(h1pbmcs.seurat.obj)
      
      #Finding the list of all post filtering commands being run on Seurat object Note: Normalization is one of the commands. Look for @ commands slot
      
      str(h1pbmcs.seurat.obj)
      
      #B3. Identifying the 10 most highly variable genes
      
      
      
  
    

    
  
    
  
    






