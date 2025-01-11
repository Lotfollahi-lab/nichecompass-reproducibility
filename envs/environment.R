if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
remotes::install_github("mojaveazure/seurat-disk")
devtools::install_github('satijalab/seurat-data')

if (!requireNamespace("SRTsim", quietly = TRUE)) {
  install.packages("SRTsim")
}

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("rhdf5")
BiocManager::install("Seurat")
BiocManager::install("scRNAseq")
BiocManager::install("zellkonverter")