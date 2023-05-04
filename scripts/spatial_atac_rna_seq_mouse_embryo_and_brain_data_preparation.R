setwd("~/workspace/projects/autotalker-reproducibility")

install.packages('Seurat')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

devtools::install_github('satijalab/seurat-data')

library(Seurat)
library(SeuratData)
library(SeuratDisk)

#####################
##### Mouse E13 #####
#####################

mouse_embryo_e13 <- readRDS(file = "datasets/srt_data/bronze/spatial_atac_rna_seq_mouse_embryo_and_brain/E13_spatial_RNA_ATAC.rds")

SaveH5Seurat(mouse_embryo_e13, filename = "datasets/srt_data/silver/spatial_atac_rna_seq_mouse_e13.h5Seurat")
Convert("datasets/srt_data/silver/spatial_atac_rna_seq_mouse_e13.h5Seurat", dest = "h5ad")