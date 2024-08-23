import os
import pathlib
import anndata as ad
import squidpy as sq
import scanpy as sc

config = {
    "n_spatially_variable_genes": 5000,
    "min_cell_gene_thresh_ratio": 0.0005,
}


def filter_spatially_variable_genes(adata, n_svg, min_cell_gene_thresh_ratio):

    # Filter genes
    min_cells = int(adata.shape[0] * min_cell_gene_thresh_ratio)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Identify spatially variable genes
    sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names)
    svg_genes = adata.uns["moranI"].index[:n_svg].tolist()
    adata.var["spatially_variable"] = adata.var_names.isin(svg_genes)
    filtered_anndata = adata[:, adata.var["spatially_variable"] == True]

    return filtered_anndata


dataset = "dataset-with-neighborhood-graph"
anndata = ad.read_h5ad(os.path.join(dataset, "dataset_with_neighborhood_graph.h5ad"))
packaged_anndata = filter_spatially_variable_genes(anndata, config["n_spatially_variable_genes"], config["min_cell_gene_thresh_ratio"])
os.makedirs("dataset-with-spatially-variable-genes", exist_ok=True)
file_name = "dataset_with_spatially_variable_genes.h5ad"
local_file = os.path.join("dataset-with-spatially-variable-genes", file_name)
packaged_anndata.write_h5ad(pathlib.Path(local_file))
