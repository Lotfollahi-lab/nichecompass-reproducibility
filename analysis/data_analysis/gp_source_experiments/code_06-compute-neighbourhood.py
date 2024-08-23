import os
import anndata as ad
import squidpy as sq
import scanpy as sc
import scipy.sparse as sp
import pathlib

config = {
   "n_neighbors": 8,
}


def load(anndata_directory):
    anndata_list = []
    for anndata in os.listdir(anndata_directory):
        adata = ad.read(os.path.join(anndata_directory, anndata))
        anndata_list.append(adata)

    return anndata_list


def compute_neighbourhood(anndata_directory, n_neighs):
    adj_key = "spatial_connectivities"
    adata_batch_list = []

    for batch in os.listdir(anndata_directory):
        print(f"Processing batch {batch}...")
        print("Loading data...")
        adata_batch = sc.read_h5ad(os.path.join(anndata_directory, batch))

        print("Computing spatial neighborhood graph...\n")
        # Compute (separate) spatial neighborhood graphs
        sq.gr.spatial_neighbors(adata_batch,
                                coord_type="generic",
                                n_neighs=n_neighs)

        # Make adjacency matrix symmetric
        adata_batch.obsp[adj_key] = (
            adata_batch.obsp[adj_key].maximum(
                adata_batch.obsp[adj_key].T))
        adata_batch_list.append(adata_batch)

    adata = ad.concat(adata_batch_list, join="inner", label="sample")

    # Combine spatial neighborhood graphs as disconnected components
    batch_connectivities = []
    len_before_batch = 0
    for i in range(len(adata_batch_list)):
        if i == 0:  # first batch
            after_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[0].shape[0],
                 (adata.shape[0] -
                  adata_batch_list[0].shape[0])))
            batch_connectivities.append(sp.hstack(
                (adata_batch_list[0].obsp[adj_key],
                 after_batch_connectivities_extension)))
        elif i == (len(adata_batch_list) - 1):  # last batch
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                 (adata.shape[0] -
                  adata_batch_list[i].shape[0])))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                 adata_batch_list[i].obsp[adj_key])))
        else:  # middle batches
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0], len_before_batch))
            after_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                 (adata.shape[0] -
                  adata_batch_list[i].shape[0] -
                  len_before_batch)))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                 adata_batch_list[i].obsp[adj_key],
                 after_batch_connectivities_extension)))
        len_before_batch += adata_batch_list[i].shape[0]
    adata.obsp[adj_key] = sp.vstack(batch_connectivities)

    return adata


raw_dataset = "mouse-organogenesis-transcriptomics"
connected_anndata = compute_neighbourhood(raw_dataset, config["n_neighbors"])
os.makedirs("dataset-with-neighborhood-graph", exist_ok=True)
file_name = "dataset_with_neighborhood_graph.h5ad"
local_file = os.path.join("dataset-with-neighborhood-graph", file_name)
connected_anndata.write_h5ad(pathlib.Path(local_file))
