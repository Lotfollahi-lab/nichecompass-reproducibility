if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("zellkonverter")
BiocManager::install("scRNAseq")

library(data.table)
library(SingleCellExperiment)
library(zellkonverter)

counts <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/counts.Rds","rb"))
metadata <- readRDS(url("https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/metadata.Rds","rb"))

sce <- SingleCellExperiment(
  list(counts=counts),
  colData=DataFrame(Area=metadata$Area,
                    celltype_mapped_refined=metadata$celltype_mapped_refined))

writeH5AD(sce, file = "seqfish_mouse_organogenesis.h5ad")