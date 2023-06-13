library(SummarizedExperiment)
library(SeuratDisk)
library(dyno)
library(tidyverse)
library(anndata)
library(Seurat)
library(loomR)

library(zellkonverter)

pbmc.loom <- as.loom(pbmc, filename = "./data/processed/adata_annotated.loom", verbose = FALSE)
# seurat_file = as.Seurat("../data/processed/adata_annotated.loom", 
#                         cells="obs_names", 
#                         features="var_names", 
#                         normalized="/matrix")

