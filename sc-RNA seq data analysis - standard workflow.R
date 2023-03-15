#script to perform standard workflow steps to analyze single cell RNA-Seq data
#data: 10 Human PBMCs
#data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728

setwd('C:/Users/ermin/Documents/Bioinformatics Shift Project/Website Portfolio Project/Bioinformatics Portfolio Website/genomics/genomics_project1/GSE150728')

#loading libraries
library(Seurat)
library(tidyverse)

#PERFORMING STANDARD WORKFLOW ON sc-RNA DATA FROM PBMCS OF A HEALTHY PATIENT LABELED 'H1'

#loading the NSCLC dataset from downloaded RDS file
nsclc.sparse.m <- readRDS('C:/Users/ermin/Documents/Bioinformatics Shift Project/Website Portfolio Project/Bioinformatics Portfolio Website/genomics/genomics_project1/GSE150728/Raw Data/GSE150728_RAW/H1/GSM4557334_HIP002_cell.counts.matrices.rds')

#running str command to see all fields in the class
str(nsclc.sparse.m)

#selecting gene expression data and storing it into variable 'h1.counts.exon'
h1.counts.exon <- nsclc.sparse.m$exon

#display first 10 rows of the h1.counts.exon matrix
head(h1.counts.exon, n=10)


