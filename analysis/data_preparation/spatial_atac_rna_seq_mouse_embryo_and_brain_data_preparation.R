setwd("~/workspace/projects/nichecompass-reproducibility")

install.packages('Seurat')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

devtools::install_github('satijalab/seurat-data')

library(Seurat)
library(SeuratData)
library(SeuratDisk)

############################
##### Mouse Embryo E13 #####
############################

mouse_e13 <- readRDS(file = "datasets/srt_data/bronze/spatial_atac_rna_seq_mouse_embryo/E13_spatial_RNA_ATAC.rds")

SaveH5Seurat(mouse_e13, filename = "datasets/srt_data/silver/spatial_atac_rna_seq_mouse_embryo_e13.h5Seurat")
Convert("datasets/srt_data/silver/spatial_atac_rna_seq_mouse_embryo_e13.h5Seurat", dest = "h5ad")

###########################
##### Mouse Brain P22 #####
###########################

mouse_p22 <- readRDS(file = "datasets/srt_data/bronze/spatial_atac_rna_seq_mouse_brain/P22mousebrain_spatial_RNA_ATAC.rds")

SaveH5Seurat(mouse_p22, filename = "datasets/srt_data/silver/spatial_atac_rna_seq_mouse_brain_p22.h5Seurat")
Convert("datasets/srt_data/silver/spatial_atac_rna_seq_mouse_brain_p22.h5Seurat", dest = "h5ad")