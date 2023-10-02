print('imports')
import os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# Define paths
srt_data_folder_path = "../../datasets/srt_data" # spatially resolved transcriptomics data
srt_data_bronze_folder_path = f"{srt_data_folder_path}/bronze"
srt_data_silver_folder_path = f"{srt_data_folder_path}/silver"
srt_data_gold_folder_path = f"{srt_data_folder_path}/gold"

# Create required directories
os.makedirs(srt_data_bronze_folder_path, exist_ok=True)
os.makedirs(srt_data_silver_folder_path, exist_ok=True)
os.makedirs(srt_data_gold_folder_path, exist_ok=True)

dataset = "nanostring_cosmx_human_nsclc"
cell_type_key = "cell_type"

batches = ["Lung5_Rep1",
           "Lung5_Rep2",
           "Lung5_Rep3",
           "Lung6",
           "Lung9_Rep1",
           "Lung9_Rep2",
           "Lung12",
           "Lung13"]

annotation_df = pd.read_csv(f"{srt_data_bronze_folder_path}/{dataset}/metadata_giotto.csv", index_col=0)

for batch_idx, batch in enumerate(batches):
    print(batch)
    gene_expr_df = pd.read_csv(f"{srt_data_bronze_folder_path}/{dataset}/{batch}/{batch}-Flat_files_and_images/{batch}_exprMat_file.csv")
    metadata_df = pd.read_csv(f"{srt_data_bronze_folder_path}/{dataset}/{batch}/{batch}-Flat_files_and_images/{batch}_metadata_file.csv")

    adata = ad.AnnData(gene_expr_df[gene_expr_df.columns.difference(["fov", "cell_ID"])].values,
                       obs=gene_expr_df[["fov", "cell_ID"]],
                       dtype="float32")
    adata.var_names = gene_expr_df.columns.difference(["fov", "cell_ID"])
    adata.obs["batch"] = batch.lower()
    
    # Add spatial coordinates from metadata
    adata.obs = pd.merge(adata.obs, metadata_df, on=["fov", "cell_ID"], how="left")
    adata.obsm["spatial"] = np.array(adata.obs[["CenterX_global_px", "CenterY_global_px"]])
    
    # Drop obs without metadata
    adata.obs.reset_index(drop=True, inplace=True)
    adata = adata[adata.obs.index.isin(adata.obs.dropna().index), :].copy()
    
    # Add cell type annotations, remove cells without annotations, and make fov unique across batches
    adata.obs["cell_ID"] = f"c_{batch_idx + 1}_" + adata.obs["fov"].astype("str") + "_" + adata.obs["cell_ID"].astype("str")
    adata.obs = pd.merge(adata.obs, annotation_df, on="cell_ID", how="left")
    adata.obs["fov"] = adata.obs["batch"] + "_" + adata.obs["fov_x"].astype(str)
    adata.obs = adata.obs[["cell_ID", "patient", "batch", "fov", "cell_type", "niche"]]
    adata.obs.index = adata.obs.index.astype(str)
    adata = adata[adata.obs.index.isin(adata.obs.dropna().index), :].copy()
    adata.obs.reset_index(drop=True, inplace=True)
    
    # Convert cell type annotations to coarser resolution
    adata.obs["cell_type_original"] = adata.obs["cell_type"].astype(str)
    adata.obs.loc[adata.obs["cell_type_original"].apply(
        lambda x: "T" in x),"cell_type"] = "NK/T cell"
    adata.obs.loc[adata.obs["cell_type_original"].apply(
        lambda x: "tumor" in x),"cell_type"] = "tumor"
    adata.obs.loc[adata.obs["cell_type_original"] == "NK","cell_type"] = "NK/T cell"
    adata.obs.loc[adata.obs["cell_type_original"].apply(
        lambda x: "DC" in x),"cell_type"] = "DC"
    adata.obs.loc[(adata.obs["cell_type_original"] == "monocyte") |
                  (adata.obs["cell_type_original"] == "macrophage") |
                  (adata.obs["cell_type"] == "DC"),"cell_type"] = "myeloid"
    
    # Remove negative probes
    adata.var.index = adata.var.index.map(str)
    adata = adata[:, ~adata.var_names.str.contains("NegPrb")].copy()
    
    # Remove low quality cells
    sc.pp.filter_cells(adata, min_counts=50)
        
    # Store gene expression in sparse format
    adata.X = sp.csr_matrix(np.array(adata.X))
    adata.layers["counts"] = adata.X.copy()
    
    # Log normalize counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    adata.write(f"{srt_data_gold_folder_path}/{dataset}_batch{batch_idx + 1}.h5ad")

# Filter genes with strong fov effects

for batch_idx, batch in enumerate(batches):
    print(batch)
    adata = sc.read_h5ad(f"{srt_data_gold_folder_path}/{dataset}_batch{batch_idx + 1}.h5ad")

    filter_genes = [
        'AKT1', 'ANXA2', 'B2M', 'C1QC', 'CALM2', 'CD3E', 'CD63', 'CD74', 'CD81',
        'CLDN4', 'COL1A1', 'COL3A1', 'COL6A1', 'COL6A2', 'COL9A2', 'CRP',
        'CTNNB1', 'DCN', 'DUSP5', 'EIF5A', 'FKBP11', 'GLUL', 'GPNMB', 'GPX1',
        'GSTP1', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB5', 'HLA-E', 'HSP90AA1',
        'HSP90AB1', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPB1', 'IFI27', 'IFITM3',
        'IGFBP5', 'IGHG1', 'IGHG2', 'IGKC', 'ITGB8', 'JUN', 'JUNB', 'KRT14',
        'KRT19', 'KRT8', 'MALAT1', 'MIF', 'MT2A', 'MZT2A', 'NDRG1', 'NEAT1',
        'RPL21', 'RPL22', 'RPL32', 'RPL34', 'RPL37', 'S100A10', 'S100A4',
        'S100A6', 'S100A9', 'SAT1', 'SERPINA1', 'SLPI', 'SOD2', 'SOX9',
        'SQSTM1', 'TACSTD2', 'TAGLN', 'THBS2', 'TIMP1', 'TM4SF1', 'TUBB',
        'TYK2', 'VEGFA']

    gene_mask = [gene not in filter_genes for gene in adata.var_names]

    adata = adata[:, gene_mask]
    adata.write(f"{srt_data_gold_folder_path}/{dataset}_filtered_batch{batch_idx + 1}.h5ad")
