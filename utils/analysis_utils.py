import collections
import colorsys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from matplotlib import gridspec



def plot_latent(adata,
                dataset_label,
                color_by,
                color_palette,
                groups,
                save_fig,
                file_path):
    """Plot features in latent space."""
    fig = sc.pl.umap(adata,
                     color=[color_by],
                     palette=color_palette,
                     groups=groups,
                     size=320000/len(adata),
                     return_fig=True)
    fig.set_size_inches(12, 8)
    plt.title(f"{dataset_label}: NicheCompass Latent "
              f"{color_by.replace('_', ' ').title()} Annotations",
              size=20,
              pad=15)

    if save_fig:
        fig.savefig(file_path,
                    bbox_inches="tight")
    else:
        fig.show()
        
        
def plot_category_in_latent_and_physical_space(
        adata,
        plot_label,
        cat_key,
        groups,
        sample_key,
        samples,
        cat_colors,
        size,
        spot_size,
        save_fig,
        file_path,
        add_outline=False):
    """Plot latent clusters in latent and physical space."""
    ncols = min(3, len(samples))
    # Create plot of cell type annotations in physical and latent space
    fig = plt.figure(figsize=(10, len(samples) + 20))
    title = fig.suptitle(t=f"{plot_label} in NicheCompass " \
                           "Latent and Physical Space",
                         y=0.96,
                         x=0.55,
                         fontsize=20)
    spec1 = gridspec.GridSpec(ncols=1,
                              nrows=1 + 1 + int(len(samples) / 3),
                              width_ratios=[1],
                              height_ratios=[3] + [2] * (1 + int(len(samples) / 3)))
    spec2 = gridspec.GridSpec(ncols=ncols,
                              nrows=2 + int(len(samples) / 3),
                              width_ratios=[1] * ncols,
                              height_ratios=[3] + [2] * (1 + int(len(samples) / 3)))
    axs = []
    axs.append(fig.add_subplot(spec1[0]))
    sc.pl.umap(adata=adata,
               color=[cat_key],
               groups=groups,
               palette=cat_colors,
               size=size,
               title=f"{plot_label} in NicheCompass Latent Space",
               add_outline=add_outline,
               ax=axs[0],
               show=False)
    for idx, sample in enumerate(samples):
        axs.append(fig.add_subplot(spec2[ncols + idx]))
        sc.pl.spatial(adata=adata[adata.obs[sample_key] == sample],
                      color=[cat_key],
                      groups=groups,                  
                      palette=cat_colors,
                      spot_size=spot_size,
                      title=f"{plot_label} in \n Physical Space \n"
                            f"(Sample: {sample})",
                      legend_loc=None,
                      ax=axs[idx+1],
                      show=False)
    
    # Create and position shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(handles,
                     labels,
                     frameon=False,
                     loc="center left",
                     bbox_to_anchor=(0.98, 0.5))

    for handle in lgd.legendHandles:
        handle.set_sizes([100])  # Adjust the size as needed
    try:
        axs[0].get_legend().remove()
    except:
        print("Could not remove axis.")

    # Adjust, save and display plot
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    if save_fig:
        fig.savefig(file_path,
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()


def plot_physical_latent_for_cell_types(adata,
                                        cell_types,
                                        sample_key,
                                        cell_type_key,
                                        cell_type_colors,
                                        figure_folder_path,
                                        save_fig=False):
    # Create plot of cell type annotations in physical and latent space
    fig = plt.figure(figsize=(12, 14))
    title = fig.suptitle(t=f"{' & '.join(cell_types).title()} in Physical and Latent Space",
                         y=0.96,
                         x=0.55,
                         fontsize=20)
    spec = gridspec.GridSpec(ncols=1,
                             nrows=2,
                             width_ratios=[1],
                             height_ratios=[1, 5])
    ax1 = fig.add_subplot(2, 3, 1)
    ax2 = fig.add_subplot(2, 3, 2)
    ax3 = fig.add_subplot(2, 3, 3)
    ax4 = fig.add_subplot(spec[1])
    sc.pl.spatial(adata=adata[adata.obs[sample_key] == "embryo1"],
                  color=[cell_type_key],
                  palette=cell_type_colors,
                  groups=cell_types,
                  size=160000/len(adata),
                  spot_size=0.03,
                  title="Physical Space Embryo 1",
                  legend_loc=None,
                  ax=ax1,
                  show=False)
    sc.pl.spatial(adata=adata[adata.obs[sample_key] == "embryo2"],
                  color=[cell_type_key],
                  palette=cell_type_colors,
                  groups=cell_types,
                  size=160000/len(adata),
                  spot_size=0.03,
                  title="Physical Space Embryo 2",
                  legend_loc=None,
                  ax=ax2,
                  show=False)
    sc.pl.spatial(adata=adata[adata.obs[sample_key] == "embryo3"],
                  color=[cell_type_key],
                  palette=cell_type_colors,
                  groups=cell_types,
                  size=160000/len(adata),
                  spot_size=0.03,
                  title="Physical Space Embryo 3",
                  legend_loc=None,
                  ax=ax3,
                  show=False)
    sc.pl.umap(adata=adata,
               color=[cell_type_key],
               palette=cell_type_colors,
               groups=cell_types,
               size=1280000/len(adata),
               title="NicheCompass Latent Space",
               ax=ax4,
               show=False)

    # Create and position shared legend
    handles, labels = ax4.get_legend_handles_labels()
    lgd = fig.legend(handles,
                     labels,
                     loc="center left",
                     bbox_to_anchor=(0.98, 0.5))
    ax4.get_legend().remove()

    # Adjust, save and display plot
    plt.subplots_adjust(wspace=0., hspace=0.85)
    if save_fig:
        fig.savefig(f"{figure_folder_path}/"
                    f"{'_'.join(cell_types).replace('/', '_').replace(' ', '_').lower()}_physical_latent_space.svg",
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()
    

def compute_latent_clusters(adata,
                            latent_resolution,
                            latent_cluster_key,
                            latent_knng_key):    

    sc.tl.leiden(adata=adata,
                 resolution=latent_resolution,
                 key_added=latent_cluster_key,
                 neighbors_key=latent_knng_key)
    
    
def compute_cell_type_latent_clusters(adata,
                                      latent_key,
                                      cell_type_latent_resolution,
                                      cell_type_res_latent_cluster_key,
                                      cell_type_latent_cluster_key,
                                      latent_knng_key,
                                      cell_type_key,
                                      cell_type,
                                      cell_type_specific_clustering=False,
                                      dummy_category="-1"):    
    if not cell_type_specific_clustering:
        # Compute latent Leiden clustering with cell-type-specific resolution
        sc.tl.leiden(adata=adata,
                     resolution=cell_type_latent_resolution,
                     key_added=cell_type_res_latent_cluster_key,
                     neighbors_key=latent_knng_key)
        
        # Filter for cell type
        cell_type_adata = adata[adata.obs[cell_type_key] == cell_type]
    else:
        # Filter for cell type
        cell_type_adata = adata[adata.obs[cell_type_key] == cell_type].copy()
        
        # Use latent representation to compute cell-type-specific nearest neighbor graph
        sc.pp.neighbors(cell_type_adata,
                        use_rep=latent_key,
                        key_added=latent_knng_key)
        
        # Compute latent Leiden clustering with cell-type-specific resolution
        sc.tl.leiden(adata=cell_type_adata,
                     resolution=cell_type_latent_resolution,
                     key_added=cell_type_res_latent_cluster_key,
                     neighbors_key=latent_knng_key)

    # Only keep latent clusters for cell type and set rest to dummy category
    adata.obs[cell_type_latent_cluster_key] = dummy_category
    adata.obs.loc[adata.obs[cell_type_key] == cell_type, cell_type_latent_cluster_key] = (
        cell_type_adata.obs[cell_type_res_latent_cluster_key])
    
    
def plot_cell_type_latent_clusters(adata,
                                   cell_type_latent_cluster_key,
                                   cell_type,
                                   groups,
                                   condition_key,
                                   latent_cluster_colors,
                                   save_fig,
                                   file_path):
    """Plot cell-type-specific latent clusters."""
    # Create plot of cell type annotations in physical and latent space
    fig = plt.figure(figsize=(12, 14))
    title = fig.suptitle(t=f"{cell_type.replace('_', ' ').title()} Latent Clusters in Latent and Physical Space",
                         y=0.96,
                         x=0.55,
                         fontsize=20)
    spec1 = gridspec.GridSpec(ncols=1,
                              nrows=2,
                              width_ratios=[1],
                              height_ratios=[3, 2])
    spec2 = gridspec.GridSpec(ncols=3,
                              nrows=2,
                              width_ratios=[1, 1, 1],
                              height_ratios=[3, 2])
    ax1 = fig.add_subplot(spec1[0])
    ax2 = fig.add_subplot(spec2[3])
    ax3 = fig.add_subplot(spec2[4])
    ax4 = fig.add_subplot(spec2[5])
    sc.pl.umap(adata=adata,
               color=[cell_type_latent_cluster_key],
               groups=groups,
               palette=latent_cluster_colors,
               size=1280000/len(adata),
               title=f"{cell_type.replace('_', ' ').title()} Latent Clusters in Latent Space",
               ax=ax1,
               show=False)
    sc.pl.spatial(adata=adata[adata.obs[condition_key] == "embryo1"],
                  color=[cell_type_latent_cluster_key],
                  groups=groups,                  
                  palette=latent_cluster_colors,
                  size=160000/len(adata),
                  spot_size=0.03,
                  title=f"{cell_type.replace('_', ' ').title()} \n Latent Clusters in \n Physical Space Embryo 1",
                  legend_loc=None,
                  ax=ax2,
                  show=False)
    sc.pl.spatial(adata=adata[adata.obs[condition_key] == "embryo2"],
                  color=[cell_type_latent_cluster_key],
                  groups=groups,
                  palette=latent_cluster_colors,
                  size=160000/len(adata),
                  spot_size=0.03,
                  title=f"{cell_type.replace('_', ' ').title()} \n Latent Clusters in \n Physical Space Embryo 2",
                  legend_loc=None,
                  ax=ax3,
                  show=False)
    sc.pl.spatial(adata=adata[adata.obs[condition_key] == "embryo3"],
                  color=[cell_type_latent_cluster_key],
                  groups=groups,
                  palette=latent_cluster_colors,
                  size=160000/len(adata),
                  spot_size=0.03,
                  title=f"{cell_type.replace('_', ' ').title()} \n Latent Clusters in \n Physical Space Embryo 3",
                  legend_loc=None,
                  ax=ax4,
                  show=False)

    # Create and position shared legend
    handles, labels = ax1.get_legend_handles_labels()
    lgd = fig.legend(handles,
                     labels,
                     loc="center left",
                     bbox_to_anchor=(0.98, 0.5))
    ax1.get_legend().remove()

    # Adjust, save and display plot
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    if save_fig:
        fig.savefig(file_path,
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()


def generate_gp_info_plots(analysis_label,
                           differential_gp_test_results_key,
                           model,
                           cell_type_key,
                           cell_type_colors,
                           latent_cluster_colors,
                           plot_category,
                           log_bayes_factor_thresh,
                           n_top_enriched_gps=10,
                           adata=None,
                           feature_spaces=["latent", "physical_embryo1", "physical_embryo2", "physical_embryo3"],
                           plot_types=["gene_categories", "top_genes"],
                           n_top_genes_per_gp=3,
                           save_figs=False,
                           figure_folder_path="",
                           spot_size=30):
    
    if adata is None:
        adata = model.adata.copy()
        
    cats = adata.uns[differential_gp_test_results_key]["category"][:n_top_enriched_gps]
    gps = adata.uns[differential_gp_test_results_key]["gene_program"][:n_top_enriched_gps]
    
    for gp in gps:
        gp_gene_importances_df = model.compute_gp_gene_importances(selected_gp=gp)
        
        if "gene_categories" in plot_types:
            pos_sign_target_genes = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] > 0) &
                (gp_gene_importances_df["gene_entity"] == "target"), "gene"].tolist()
            pos_sign_source_genes = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] > 0) &
                (gp_gene_importances_df["gene_entity"] == "source"), "gene"].tolist()
            neg_sign_target_genes = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] < 0) &
                (gp_gene_importances_df["gene_entity"] == "target"), "gene"].tolist()
            neg_sign_source_genes = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] < 0) &
                (gp_gene_importances_df["gene_entity"] == "source"), "gene"].tolist()

            pos_sign_target_gene_importances = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] > 0) &
                (gp_gene_importances_df["gene_entity"] == "target"), "gene_importance"].values.reshape(1, -1)
            pos_sign_source_gene_importances = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] > 0) &
                (gp_gene_importances_df["gene_entity"] == "source"), "gene_importance"].values.reshape(1, -1)
            neg_sign_target_gene_importances = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] < 0) &
                (gp_gene_importances_df["gene_entity"] == "target"), "gene_importance"].values.reshape(1, -1)
            neg_sign_source_gene_importances = gp_gene_importances_df.loc[
                (gp_gene_importances_df["gene_weight_sign_corrected"] < 0) &
                (gp_gene_importances_df["gene_entity"] == "source"), "gene_importance"].values.reshape(1, -1)

            pos_sign_target_gene_expr = adata[:, pos_sign_target_genes].X.toarray()
            pos_sign_source_gene_expr = adata[:, pos_sign_source_genes].X.toarray()
            neg_sign_target_gene_expr = adata[:, neg_sign_target_genes].X.toarray()
            neg_sign_source_gene_expr = adata[:, neg_sign_source_genes].X.toarray()

            adata.obs[f"{gp}_pos_sign_target_gene_weighted_mean_gene_expr"] = (
                np.mean(pos_sign_target_gene_expr * pos_sign_target_gene_importances, axis=1))
            adata.obs[f"{gp}_pos_sign_source_gene_weighted_mean_gene_expr"] = (
                np.mean(pos_sign_source_gene_expr * pos_sign_source_gene_importances, axis=1))
            adata.obs[f"{gp}_neg_sign_target_gene_weighted_mean_gene_expr"] = (
                np.mean(neg_sign_target_gene_expr * neg_sign_target_gene_importances, axis=1))
            adata.obs[f"{gp}_neg_sign_source_gene_weighted_mean_gene_expr"] = (
                np.mean(neg_sign_source_gene_expr * neg_sign_source_gene_importances, axis=1))

            adata.uns[f"{gp}_gene_category_importances"] = np.array([pos_sign_target_gene_importances.sum(),
                                                                     pos_sign_source_gene_importances.sum(),
                                                                     neg_sign_target_gene_importances.sum(),
                                                                     neg_sign_source_gene_importances.sum()])
        
        if "top_genes" in plot_types:
            gp_source_genes_gene_importances_df = gp_gene_importances_df[
                gp_gene_importances_df["gene_entity"] == "source"]
            
            gp_target_genes_gene_importances_df = gp_gene_importances_df[
                gp_gene_importances_df["gene_entity"] == "target"]
            
            adata.uns["n_top_source_genes"] = n_top_genes_per_gp
            adata.uns[f"{gp}_source_genes_top_genes"] = gp_source_genes_gene_importances_df["gene"][:n_top_genes_per_gp]
            adata.uns[f"{gp}_source_genes_top_gene_importances"] = gp_source_genes_gene_importances_df["gene_importance"][:n_top_genes_per_gp]
            adata.uns[f"{gp}_source_genes_top_gene_signs"] = np.where(gp_source_genes_gene_importances_df["gene_weight_sign_corrected"] > 0, "+", "-")
            adata.uns["n_top_target_genes"] = n_top_genes_per_gp
            adata.uns[f"{gp}_target_genes_top_genes"] = gp_target_genes_gene_importances_df["gene"][:n_top_genes_per_gp]
            adata.uns[f"{gp}_target_genes_top_gene_importances"] = gp_target_genes_gene_importances_df["gene_importance"][:n_top_genes_per_gp]
            adata.uns[f"{gp}_target_genes_top_gene_signs"] = np.where(gp_target_genes_gene_importances_df["gene_weight_sign_corrected"] > 0, "+", "-")
            
            #adata.uns["n_top_genes"] = n_top_genes_per_gp
            #adata.uns[f"{gp}_top_genes"] = gp_gene_importances_df["gene"][:n_top_genes_per_gp]
            #adata.uns[f"{gp}_top_gene_importances"] = gp_gene_importances_df["gene_importance"][:n_top_genes_per_gp]
            #adata.uns[f"{gp}_top_gene_signs"] = np.where(gp_gene_importances_df["gene_weight_sign_corrected"] > 0, "+", "-")
            #adata.uns[f"{gp}_top_gene_entities"] = gp_gene_importances_df["gene_entity"]
        
    for feature_space in feature_spaces:
        for plot_type in plot_types:
            plot_gp_info_plots(adata=adata,
                               cell_type_key=cell_type_key,
                               cell_type_colors=cell_type_colors,
                               latent_cluster_colors=latent_cluster_colors,
                               cats=cats,
                               gps=gps,
                               plot_type=plot_type,
                               plot_category=plot_category,
                               feature_space=feature_space,
                               spot_size=spot_size,
                               suptitle=f"{analysis_label.replace('_', ' ').title()} Top {n_top_enriched_gps} Enriched GPs: "
                                        f"GP Scores and {'Weighted Mean ' if plot_type == 'gene_categories' else ''}"
                                        f"Gene Expression of {plot_type.replace('_', ' ').title()} in {feature_space.replace('_', ' ').title()} Feature Space",
                               cat_title=f"Enriched GP Category in \n {plot_category.replace('_', ' ').title()}",
                               save_fig=save_figs,
                               figure_folder_path=figure_folder_path,
                               fig_name=f"{analysis_label}_log_bayes_factor_{log_bayes_factor_thresh}_enriched_gps_gp_scores_" \
                                        f"{'weighted_mean_gene_expr' if plot_type == 'gene_categories' else 'top_genes_gene_expr'}_" \
                                        f"{feature_space}_space")
            
def plot_gp_info_plots(adata,
                       cell_type_key,
                       cell_type_colors,
                       latent_cluster_colors,
                       cats,
                       gps,
                       plot_type,
                       plot_category,
                       feature_space,
                       spot_size,
                       suptitle,
                       cat_title,
                       save_fig,
                       figure_folder_path,
                       fig_name):
    if plot_category == cell_type_key:
        palette = cell_type_colors
    else:
        palette = latent_cluster_colors 
    # Plot selected gene program latent scores
    if plot_type == "gene_categories":
        ncols = 6
        fig_width = 36
        wspace = 0.155
    elif plot_type == "top_genes":
        ncols = 2 + adata.uns["n_top_genes"]
        fig_width = 12 + (6 * adata.uns["n_top_genes"])
        wspace = 0.3
    fig, axs = plt.subplots(nrows=len(gps), ncols=ncols, figsize=(fig_width, 6*len(gps)))
    if axs.ndim == 1:
        axs = axs.reshape(1, -1)

    title = fig.suptitle(t=suptitle,
                         x=0.55,
                         y=(1.1 if len(gps) == 1 else 0.93),
                         fontsize=20)
    for i, gp in enumerate(gps):
        if feature_space == "latent":
            sc.pl.umap(adata,
                       color=plot_category,
                       palette=palette,
                       groups=cats[i],
                       ax=axs[i, 0],
                       title=cat_title,
                       legend_loc="on data",
                       na_in_legend=False,
                       show=False)
            sc.pl.umap(adata,
                       color=gps[i],
                       color_map="RdBu",
                       ax=axs[i, 1],
                       title=f"{gp[:gp.index('_')]}\n{gp[gp.index('_') + 1: gp.rindex('_')].replace('_', ' ')}\n{gp[gps[i].rindex('_') + 1:]} score",
                       show=False)
        elif "physical" in feature_space:
            sc.pl.spatial(adata=adata[adata.obs["batch"] == feature_space.split("_")[1]],
                          color=plot_category,
                          palette=palette,
                          groups=cats[i],
                          ax=axs[i, 0],
                          spot_size=spot_size,
                          title=cat_title,
                          legend_loc="on data",
                          na_in_legend=False,
                          show=False)
            sc.pl.spatial(adata=adata[adata.obs["batch"] == feature_space.split("_")[1]],
                          color=gps[i],
                          color_map="RdBu",
                          spot_size=spot_size,
                          title=f"{gps[i].split('_', 1)[0]}\n{gps[i].split('_', 1)[1]}",
                          legend_loc=None,
                          ax=axs[i, 1],
                          show=False) 
        axs[i, 0].xaxis.label.set_visible(False)
        axs[i, 0].yaxis.label.set_visible(False)
        axs[i, 1].xaxis.label.set_visible(False)
        axs[i, 1].yaxis.label.set_visible(False)
        if plot_type == "gene_categories":
            for j, gene_category in enumerate(["pos_sign_target_gene",
                                               "pos_sign_source_gene",
                                               "neg_sign_target_gene",
                                               "neg_sign_source_gene"]):
                if not adata.obs[f"{gp}_{gene_category}_weighted_mean_gene_expr"].isna().any():
                    if feature_space == "latent":
                        sc.pl.umap(adata,
                                   color=f"{gp}_{gene_category}_weighted_mean_gene_expr",
                                   color_map=("Blues" if "pos_sign" in gene_category else "Reds"),
                                   ax=axs[i, j+2],
                                   legend_loc="on data",
                                   na_in_legend=False,
                                   title=f"Weighted mean gene expression \n {gene_category.replace('_', ' ')} ({adata.uns[f'{gp}_gene_category_importances'][j]:.2f})",
                                   show=False)
                    elif "physical" in feature_space:
                        sc.pl.spatial(adata=adata[adata.obs["sample"] == feature_space.split("_")[1]],
                                      color=f"{gp}_{gene_category}_weighted_mean_gene_expr",
                                      color_map=("Blues" if "pos_sign" in gene_category else "Reds"),
                                      ax=axs[i, 2+j],
                                      legend_loc="on data",
                                      na_in_legend=False,
                                      groups=cats[i],
                                      spot_size=spot_size,
                                      title=f"Weighted mean gene expression \n {gene_category.replace('_', ' ')} ({adata.uns[f'{gp}_gene_category_importances'][j]:.2f})",
                                      show=False)                        
                    axs[i, j+2].xaxis.label.set_visible(False)
                    axs[i, j+2].yaxis.label.set_visible(False)
                else:
                    axs[i, j+2].set_visible(False)
        elif plot_type == "top_genes":
            for j in range(len(adata.uns[f"{gp}_top_genes"])):
                if feature_space == "latent":
                    sc.pl.umap(adata,
                               color=adata.uns[f"{gp}_top_genes"][j],
                               color_map=("Blues" if adata.uns[f"{gp}_top_gene_signs"][j] == "+" else "Reds"),
                               ax=axs[i, 2+j],
                               legend_loc="on data",
                               na_in_legend=False,
                               title=f"{adata.uns[f'{gp}_top_genes'][j]}: "
                                     f"{adata.uns[f'{gp}_top_gene_importances'][j]:.2f} "
                                     f"({adata.uns[f'{gp}_top_gene_entities'][j][0]}; "
                                     f"{adata.uns[f'{gp}_top_gene_signs'][j]})",
                               show=False)
                elif "physical" in feature_space:
                    sc.pl.spatial(adata=adata[adata.obs["batch"] == feature_space.split("_")[1]],
                                  color=adata.uns[f"{gp}_top_genes"][j],
                                  color_map=("Blues" if adata.uns[f"{gp}_top_gene_signs"][j] == "+" else "Reds"),
                                  legend_loc="on data",
                                  na_in_legend=False,
                                  ax=axs[i, 2+j],
                                  # groups=cats[i],
                                  spot_size=spot_size,
                                  title=f"{adata.uns[f'{gp}_top_genes'][j]}: "
                                        f"{adata.uns[f'{gp}_top_gene_importances'][j]:.2f} "
                                        f"({adata.uns[f'{gp}_top_gene_entities'][j][0]}; "
                                        f"{adata.uns[f'{gp}_top_gene_signs'][j]})",
                                  show=False)
                axs[i, 2+j].xaxis.label.set_visible(False)
                axs[i, 2+j].yaxis.label.set_visible(False)
            for k in range(len(adata.uns[f"{gp}_top_genes"]), ncols - 2):
                axs[i, 2+k].set_visible(False)

    # Save and display plot
    plt.subplots_adjust(wspace=wspace, hspace=0.275)
    if save_fig:
        fig.savefig(f"{figure_folder_path}/{fig_name}.svg",
                    bbox_extra_artists=(title,),
                    bbox_inches="tight")
    plt.show()
    
    
def add_sub_cell_type(row,
                      cell_type_key,
                      cell_type,
                      cell_type_latent_cluster_key):  
    if row[cell_type_key] == cell_type:
        return row[cell_type_key] + " cluster " + row[cell_type_latent_cluster_key]
    return row[cell_type_key]


def add_cell_type_latent_cluster_emphasis(row,
                                          cell_type_latent_cluster_key,
                                          selected_latent_clusters,
                                          comparison_latent_clusters):
    """Add cell type latent cluster emphasis column for plotting."""
    for selected_latent_cluster in selected_latent_clusters:
        if row[cell_type_latent_cluster_key] == selected_latent_cluster:
            return row[cell_type_latent_cluster_key]
    for comparison_latent_cluster in comparison_latent_clusters:
        if row[cell_type_latent_cluster_key] == comparison_latent_cluster:
            return "-1"
    return np.nan


def store_top_gps_summary(model,
                          top_gps,
                          file_path):
    gp_summary_df = model.get_gp_summary()
    
    # Retrieve summary information for top gene programs
    top_gps_summary_df = gp_summary_df[gp_summary_df["gp_name"].isin(top_gps)][[
        "gp_name",
        "n_source_genes",
        "n_non_zero_source_genes",
        "n_target_genes",
        "n_non_zero_target_genes",
        "gp_source_genes",
        "gp_target_genes",
        "gp_source_genes_weights_sign_corrected",
        "gp_target_genes_weights_sign_corrected",
        "gp_source_genes_importances",
        "gp_target_genes_importances"]]

    # Write to disk
    top_gps_summary_df.to_csv(f"{file_path}")
    return top_gps_summary_df

def create_differential_gps_dotplot(adata,
                                    groupby_key,
                                    title,
                                    save_fig,
                                    file_path):
    
    # Retrieve differential gene programs
    differential_gps = [col for col in model.adata.obs.columns if col.endswith("_GP")]
    
    # Create grouped dotplot
    fig = sc.pl.dotplot(adata,
                        differential_gps,
                        groupby=groupby_key,
                        dendrogram=True, 
                        title=title,
                        swap_axes=True,
                        return_fig=True)
    
    # Save and display plot
    if save_fig:
        fig.savefig(file_path)
    fig.show()

    
def get_distinct_colors(n):
    """
    https://www.quora.com/How-do-I-generate-n-visually-distinct-RGB-colours-in-Python/answer/Karthik-Kumar-Viswanathan
    """
    hue_partition = 1 / (n + 1)
    colors = [
        colorsys.hsv_to_rgb(hue_partition * value, 1.0, 1.0) for value in range(0, n)
    ]
    return colors[::2] + colors[1::2]


def text_width(fig, ax, text, fontsize):
    text = ax.text(-100, 0, text, fontsize=fontsize)
    text_bb = text.get_window_extent(renderer=fig.canvas.get_renderer())
    text_bb = text_bb.transformed(fig.dpi_scale_trans.inverted())
    width = text_bb.width
    text.remove()
    return width
    
    
class Sankey:
    def __init__(
        self,
        x,
        y,
        colorside,
        plot_width=8,
        plot_height=8,
        gap=0.12,
        alpha=0.3,
        fontsize="small",
        left_order=None,
        mapping=None,
        colors=None,
        #                  colorside=None,
        tag=None,
        title=None,
        title_left=None,
        title_right=None,
        ax=None,
    ):
        self.X = x
        self.Y = y
        if ax:
            self.plot_width = ax.get_position().width * ax.figure.get_size_inches()[0]
            self.plot_height = ax.get_position().height * ax.figure.get_size_inches()[1]
        else:
            self.plot_width = plot_width
            self.plot_height = plot_height
        self.gap = gap
        self.alpha = alpha
        self.colors = colors
        self.colorside = colorside
        self.fontsize = fontsize
        self.tag = tag
        self.map = mapping is not None
        self.mapping = mapping
        self.mapping_colors = {
            "increase": "#1f721c",
            "decrease": "#ddc90f",
            "mistake": "#dd1616",
            "correct": "#dddddd",
            "novel": "#59a8d6",
        }
        self.title = title
        self.title_left = title_left
        self.title_right = title_right

        self.need_title = any(
            map(lambda x: x is not None, (title, title_left, title_right))
        )
        if self.need_title:
            self.plot_height -= 0.5

        self.init_figure(ax)

        self.flows = collections.Counter(zip(x, y))
        self.init_nodes(left_order)

        self.init_widths()
        # inches per 1 item in x and y
        self.resolution = (plot_height - gap * (len(self.left_nodes) - 1)) / len(x)
        if self.colors == None:
            if colorside == "left":
                self.colors = {
                    name: colour
                    for name, colour in zip(
                        self.left_nodes.keys(),
                        get_distinct_colors(len(self.left_nodes)),
                    )
                }
            elif colorside == "right":
                self.colors = {
                    name: colour
                    for name, colour in zip(
                        self.right_nodes.keys(),
                        get_distinct_colors(len(self.right_nodes)),
                    )
                }
            else:
                raise ValueError(
                    "colorside argument should be set either to 'left' or 'right'. Exiting."
                )

        self.init_offsets()

    def init_figure(self, ax):
        if ax is None:
            self.fig = plt.figure()
            self.ax = plt.Axes(self.fig, [0, 0, 1, 1])
            self.fig.add_axes(self.ax)
        self.fig = ax.figure
        self.ax = ax

    def init_nodes(self, left_order):
        left_nodes = {}
        right_nodes = {}
        left_offset = 0
        for (left, right), flow in self.flows.items():
            if left in left_nodes:
                left_nodes[left] += flow
            else:
                left_nodes[left] = flow
            if right in right_nodes:
                node = right_nodes[right]
                node[0] += flow
                if flow > node[2]:
                    node[1] = left
                    node[2] = flow
            else:
                right_nodes[right] = [flow, left, flow]

        self.left_nodes = collections.OrderedDict()
        self.left_nodes_idx = {}
        if left_order is None:
            key = lambda pair: -pair[1]
        else:
            left_order = list(left_order)
            key = lambda pair: left_order.index(pair[0])

        for name, flow in sorted(left_nodes.items(), key=key):
            self.left_nodes[name] = flow
            self.left_nodes_idx[name] = len(self.left_nodes_idx)

        left_names = list(self.left_nodes.keys())
        self.right_nodes = collections.OrderedDict()
        self.right_nodes_idx = {}
        for name, node in sorted(
            right_nodes.items(),
            key=lambda pair: (left_names.index(pair[1][1]), -pair[1][2]),
        ):
            self.right_nodes[name] = node[0]
            self.right_nodes_idx[name] = len(self.right_nodes_idx)

    def init_widths(self):
        self.left_width = max(
            (
                text_width(self.fig, self.ax, node, self.fontsize)
                for node in self.left_nodes
            )
        )
        if self.title_left:
            self.left_width = max(
                self.left_width,
                text_width(self.fig, self.ax, self.title_left, self.fontsize) / 2,
            )
        self.right_width = max(
            (
                text_width(self.fig, self.ax, node, self.fontsize)
                for node in self.right_nodes
            )
        )
        if self.title_right:
            self.right_width = max(
                self.right_width,
                text_width(self.fig, self.ax, self.title_right, self.fontsize) / 2,
            )

        self.right_stop = self.plot_width - self.left_width - self.right_width
        self.middle1_stop = self.right_stop * 9 / 20
        self.middle2_stop = self.right_stop * 11 / 20

    def init_offsets(self):
        self.offsets_l = {}
        self.offsets_r = {}

        offset = 0
        for name, flow in self.left_nodes.items():
            self.offsets_l[name] = offset
            offset += flow * self.resolution + self.gap

        offset = 0
        for name, flow in self.right_nodes.items():
            self.offsets_r[name] = offset
            offset += flow * self.resolution + self.gap

    def draw_flow(self, left, right, flow, node_offsets_l, node_offsets_r, colorside):
        P = matplotlib.path.Path

        flow *= self.resolution
        left_y = self.offsets_l[left] + node_offsets_l[left]
        right_y = self.offsets_r[right] + node_offsets_r[right]
        if self.need_title:
            left_y += 0.5
            right_y += 0.5
        node_offsets_l[left] += flow
        node_offsets_r[right] += flow
        if colorside == "left":
            color = self.colors[left]
        elif colorside == "right":
            color = self.colors[right]
        if self.mapping is not None:
            color = self.mapping_colors[self.mapping.category(left, right)]

        path_data = [
            (P.MOVETO, (0, -left_y)),
            (P.LINETO, (0, -left_y - flow)),
            (P.CURVE4, (self.middle1_stop, -left_y - flow)),
            (P.CURVE4, (self.middle2_stop, -right_y - flow)),
            (P.CURVE4, (self.right_stop, -right_y - flow)),
            (P.LINETO, (self.right_stop, -right_y)),
            (P.CURVE4, (self.middle2_stop, -right_y)),
            (P.CURVE4, (self.middle1_stop, -left_y)),
            (P.CURVE4, (0, -left_y)),
            (P.CLOSEPOLY, (0, -left_y)),
        ]
        codes, verts = zip(*path_data)
        path = P(verts, codes)
        patch = matplotlib.patches.PathPatch(
            path,
            facecolor=color,
            alpha=0.9 if flow < 0.02 else self.alpha,
            edgecolor="none",
        )
        self.ax.add_patch(patch)

    def draw_label(self, label, is_left):
        nodes = self.left_nodes if is_left else self.right_nodes
        offsets = self.offsets_l if is_left else self.offsets_r
        y = offsets[label] + nodes[label] * self.resolution / 2
        if self.need_title:
            y += 0.5

        self.ax.text(
            -0.1 if is_left else self.right_stop + 0.1,
            -y,
            label,
            horizontalalignment="right" if is_left else "left",
            verticalalignment="center",
            fontsize=self.fontsize,
        )

    def draw_titles(self):
        if self.title:
            self.ax.text(
                self.right_stop / 2,
                -0.25,
                self.title,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=self.fontsize,
                fontweight="bold",
            )
        if self.title_left:
            self.ax.text(
                -0.1,
                -0.25,
                self.title_left,
                horizontalalignment="right",
                verticalalignment="center",
                fontsize=self.fontsize,
            )
        if self.title_right:
            self.ax.text(
                self.right_stop + 0.1,
                -0.25,
                self.title_right,
                horizontalalignment="left",
                verticalalignment="center",
                fontsize=self.fontsize,
            )

    def draw(self, colorside):
        node_offsets_l = collections.Counter()
        node_offsets_r = collections.Counter()

        for (left, right), flow in sorted(
            self.flows.items(),
            key=lambda pair: (
                self.left_nodes_idx[pair[0][0]],
                self.right_nodes_idx[pair[0][1]],
            ),
        ):
            self.draw_flow(left, right, flow, node_offsets_l, node_offsets_r, colorside)

        for name in self.left_nodes:
            self.draw_label(name, True)
        for name in self.right_nodes:
            self.draw_label(name, False)
        self.draw_titles()

        self.ax.axis("equal")
        self.ax.set_xlim(
            -self.left_width - self.gap, self.right_stop + self.gap + self.right_width
        )
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        for k in self.ax.spines.keys():
            self.ax.spines[k].set_visible(False)
        # plt.axis('off')
        # self.fig.set_figheight(self.plot_height)
        # self.fig.set_figwidth(self.plot_width)
        if self.tag:
            text_ax = self.fig.add_axes((0.02, 0.95, 0.05, 0.05), frame_on=False)
            text_ax.set_axis_off()
            plt.text(
                0, 0, self.tag, fontsize=30, transform=text_ax.transAxes
            )
        # plt.tight_layout()


def sankey(x, y, colorside="left", **kwargs):
    diag = Sankey(x, y, colorside, **kwargs)
    diag.draw(colorside)
    return diag.fig