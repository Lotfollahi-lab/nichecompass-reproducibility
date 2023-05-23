setwd("~/workspace/projects/nichecompass-reproducibility")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("zellkonverter")
BiocManager::install("scRNAseq")
BiocManager::install("rhdf5")

library(data.table)
library(SingleCellExperiment)
library(zellkonverter)
library(rhdf5)

metadata <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/metadata.Rds","rb"))

########################
##### seqFISH data #####
########################

counts <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/counts.Rds", "rb"))
logcounts <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/exprs.Rds", "rb"))
sce <- SingleCellExperiment(
  list(logcounts=logcounts, counts=counts),
  colData=DataFrame(Area=metadata$Area,
                    celltype_mapped_refined=metadata$celltype_mapped_refined,
                    x=metadata$x_global_affine,
                    y=metadata$y_global_affine))
writeH5AD(sce, file = "datasets/srt_data/silver/seqfish_mouse_organogenesis.h5ad")

########################
##### Imputed data #####
########################

download.file("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/imputed.h5", "imputed.h5")
imputed_row_names <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/imputed_row_names.Rds", "rb"))
imputed_col_names <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/imputed_column_names.Rds", "rb"))
logcounts_imputed <- as(h5read("datasets/srt_data/bronze/seqfish_mouse_organogenesis/imputed.h5", "/logcounts"), "dgCMatrix")
sce <- SingleCellExperiment(
  list(logcounts=logcounts_imputed),
  colData=DataFrame(Area=metadata[imputed_col_names, ]$Area,
                    celltype_mapped_refined=metadata[imputed_col_names, ]$celltype_mapped_refined,
                    x=metadata[imputed_col_names, ]$x_global_affine,
                    y=metadata[imputed_col_names, ]$y_global_affine))
rownames(sce) <- imputed_row_names
colnames(sce) <- imputed_col_names
writeH5AD(sce, file = "datasets/srt_data/silver/seqfish_mouse_organogenesis_imputed.h5ad")