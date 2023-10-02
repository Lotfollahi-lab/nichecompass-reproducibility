import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

parser = argparse.ArgumentParser()
parser.add_argument('-cvCat', '--categorical_covariate', action='append', default=None)
args = parser.parse_args()

# File name
covariates = args.categorical_covariate
covariates = '-'.join(covariates)

# Create anndata
adata = sc.read_h5ad("datasets/srt_data/gold/nanostring_cosmx_human_nsclc_all_raw.h5ad")

# Train model
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.write("datasets/srt_data/gold/nanostring_cosmx_human_nsclc_all_raw_scVI.h5ad")

