setwd("~/workspace/projects/nichecompass-reproducibility")

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

seurat <- readRDS("datasets/st_data/bronze/nanostring_cosmx_human_liver/LiverDataReleaseSeurat_newUMAP.RDS")
sce <- as.SingleCellExperiment(seurat)
writeH5AD(sce, file = "datasets/st_data/silver/nanostring_cosmx_human_liver.h5ad")