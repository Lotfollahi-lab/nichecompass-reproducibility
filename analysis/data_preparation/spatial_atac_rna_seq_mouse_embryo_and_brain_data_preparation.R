setwd("~/workspace/projects/nichecompass-reproducibility")

library(Seurat)
library(SeuratData)
library(SeuratDisk)

############################
##### Mouse Embryo E13 #####
############################

options(timeout=1200)
mouse_e13 <- readRDS(url("https://cells.ucsc.edu/brain-spatial-omics/e13/E13_spatial_RNA_ATAC.rds"))

SaveH5Seurat(mouse_e13, filename = "datasets/st_data/silver/spatial_atac_rna_seq_mouse_embryo_e13.h5Seurat")
Convert("datasets/st_data/silver/spatial_atac_rna_seq_mouse_embryo_e13.h5Seurat", dest = "h5ad")

###########################
##### Mouse Brain P22 #####
###########################

mouse_p22 <- readRDS(url("https://cells.ucsc.edu/brain-spatial-omics/p22-atac/P22mousebrain_spatial_RNA_ATAC.rds"))

SaveH5Seurat(mouse_p22, filename = "datasets/st_data/silver/spatial_atac_rna_seq_mouse_brain_p22.h5Seurat")
Convert("datasets/st_data/silver/spatial_atac_rna_seq_mouse_brain_p22.h5Seurat", dest = "h5ad")