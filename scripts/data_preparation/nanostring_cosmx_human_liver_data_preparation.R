setwd("~/workspace/projects/autotalker-reproducibility")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("zellkonverter")

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

seurat <- readRDS("datasets/srt_data/bronze/nanostring_cosmx_human_liver/LiverDataReleaseSeurat_newUMAP.RDS")
sce <- as.SingleCellExperiment(seurat)
writeH5AD(sce, file = "datasets/srt_data/silver/nanostring_cosmx_human_liver.h5ad")

seurat <- readRDS("/lustre/groups/talaveralopez/INBOX/annotation/bronze/nanostring_cosmx_human_liver/LiverDataReleaseSeurat_newUMAP.RDS")
sce <- as.SingleCellExperiment(seurat)
writeH5AD(sce, file = "/lustre/groups/talaveralopez/INBOX/annotation/silver/nanostring_cosmx_human_liver.h5ad")