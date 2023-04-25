import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import gridspec

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
               title="Autotalker Latent Space",
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
                    f"{'_'.join(cell_types).replace('/', '_').replace(' ', '_').lower()}_physical_latent_space.png",
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()
    
    
def plot_latent_physical_for_cell_type_latent_clusters(adata,
                                                       cell_type,
                                                       groups,
                                                       save_fig=False):
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
               size=1280000/len(model.adata),
               title=f"{cell_type.replace('_', ' ').title()} Latent Clusters in Latent Space",
               ax=ax1,
               show=False)
    sc.pl.spatial(adata=adata[model.adata.obs[sample_key] == "embryo1"],
                  color=[cell_type_latent_cluster_key],
                  groups=groups,                  
                  palette=latent_cluster_colors,
                  size=160000/len(model.adata),
                  spot_size=0.03,
                  title=f"{cell_type.replace('_', ' ').title()} \n Latent Clusters in \n Physical Space Embryo 1",
                  legend_loc=None,
                  ax=ax2,
                  show=False)
    sc.pl.spatial(adata=adata[model.adata.obs[sample_key] == "embryo2"],
                  color=[cell_type_latent_cluster_key],
                  groups=groups,
                  palette=latent_cluster_colors,
                  size=160000/len(model.adata),
                  spot_size=0.03,
                  title=f"{cell_type.replace('_', ' ').title()} \n Latent Clusters in \n Physical Space Embryo 2",
                  legend_loc=None,
                  ax=ax3,
                  show=False)
    sc.pl.spatial(adata=adata[model.adata.obs[sample_key] == "embryo3"],
                  color=[cell_type_latent_cluster_key],
                  groups=groups,
                  palette=latent_cluster_colors,
                  size=160000/len(model.adata),
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
        fig.savefig(f"{figure_folder_path}/{cell_type.replace('/', '_').replace(' ', '_').lower()}_latent_clusters_physical_latent_space.png",
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()
    
    
def get_differential_analysis_results(analysis_label,
                                      model,
                                      adata,
                                      cat_key,
                                      selected_cats,
                                      differential_gp_scores_key,
                                      plot_category,
                                      plot_group,
                                      random_seed,
                                      comparison_cats=None,
                                      selected_gps=None,
                                      n_top_up_gps=3,
                                      n_top_down_gps=3,
                                      feature_spaces=["latent"], # "physical_embryo1", "physical_embryo2", "physical_embryo3"
                                      save_figs=False):
    # Compute gene program enrichments and retrieve top up- and downregulated gene programs
    top_unique_gps = model.compute_differential_gp_scores(cat_key=cat_key,
                                                          adata=adata,
                                                          selected_gps=selected_gps,
                                                          selected_cats=selected_cats,
                                                          gp_scores_weight_normalization=False,
                                                          comparison_cats=comparison_cats,
                                                          n_sample=10000,
                                                          key_added=differential_gp_scores_key,
                                                          n_top_up_gps_retrieved=n_top_up_gps,
                                                          n_top_down_gps_retrieved=n_top_down_gps,
                                                          seed=random_seed)

    # Display top upregulated gene programs
    top_up_gp_df = model.adata.uns[differential_gp_scores_key][:n_top_up_gps]
    display(top_up_gp_df)

    # Display top downregulated gene programs
    top_down_gp_df = model.adata.uns[differential_gp_scores_key][-n_top_down_gps:][::-1]
    display(top_down_gp_df)

    # Filter adata for GP dotplot
    if comparison_cats != "rest":
        adata_subset = model.adata[model.adata.obs[cat_key].isin(selected_cats + comparison_cats)]
    else:
        adata_subset = model.adata
    
    # Create GP dotplot
    sc.tl.dendrogram(adata_subset, groupby=cat_key)
    fig = sc.pl.dotplot(adata_subset,
                        top_unique_gps,
                        groupby=cat_key,
                        dendrogram=True, 
                        title=f"{analysis_label.replace('_', ' ').title()} Differential GP Scores",
                        swap_axes=True,
                        return_fig=True)
    # Save and display plot
    if save_figs:
        fig.savefig(f"{figure_folder_path}/{analysis_label}_differential_gp_scores.png")
    fig.show()

    # Inspect top up- and downregulated gene programs
    display(gp_summary_df[gp_summary_df["gp_name"].isin(top_unique_gps)])

    top_cats = top_up_gp_df["category"].append(top_down_gp_df["category"]).to_list()
    top_gps = top_up_gp_df["gene_program"].append(top_down_gp_df["gene_program"]).to_list()

    top_cats = top_up_gp_df["category"].append(top_down_gp_df["category"]).to_list()
    top_gps = top_up_gp_df["gene_program"].append(top_down_gp_df["gene_program"]).to_list()

    get_cat_gp_score_gene_expr_summary(analysis_label=analysis_label,
                                       model=model,
                                       cats=top_cats,
                                       gps=top_gps,
                                       plot_category=plot_category,
                                       plot_group=plot_group,
                                       adata=adata,
                                       feature_spaces=feature_spaces,
                                       plot_types=["gene_categories", "individual_genes"])
    
    return top_unique_gps


def get_cat_gp_score_gene_expr_summary(analysis_label,
                                       model,
                                       cats,
                                       gps,
                                       plot_category,
                                       plot_group,
                                       adata=None,
                                       feature_spaces=["latent", "physical_embryo1", "physical_embryo2", "physical_embryo3"],
                                       plot_types=["gene_categories", "individual_genes"],
                                       n_top_genes_per_gp=5,
                                       save_figs=False):
    
    if adata is None:
        adata = model.adata.copy()
    
    for gp in gps:
        gp_gene_importances_df = model.compute_gp_gene_importances(selected_gp=gp)
        
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
        
        adata.uns["n_top_genes"] = n_top_genes_per_gp
        adata.uns[f"{gp}_top_genes"] = gp_gene_importances_df["gene"][:n_top_genes_per_gp]
        adata.uns[f"{gp}_top_gene_importances"] = gp_gene_importances_df["gene_importance"][:n_top_genes_per_gp]
        adata.uns[f"{gp}_top_gene_signs"] = np.where(gp_gene_importances_df["gene_weight_sign_corrected"] > 0, "+", "-")
        adata.uns[f"{gp}_top_gene_entities"] = gp_gene_importances_df["gene_entity"]
        
    for feature_space in feature_spaces:
        for plot_type in plot_types:
            plot_cat_gp_score_gene_expr_summary(adata=adata,
                                                cats=cats,
                                                gps=gps,
                                                plot_type=plot_type,
                                                plot_category=plot_category,
                                                plot_group=plot_group,
                                                feature_space=feature_space,
                                                suptitle=f"{analysis_label.replace('_', ' ').title()} Differential GPs: "
                                                         f"GP Scores and {'Weighted Mean ' if plot_type == 'gene_categories' else ''}"
                                                         f"Gene Expression of {plot_type.replace('_', ' ').title()} in {feature_space.replace('_', ' ').title()} Space",
                                                cat_title=f"Differential GP \n {plot_category.replace('_', ' ').title()}",
                                                save_fig=save_figs,
                                                fig_name=f"{analysis_label}_gp_scores_{'weighted_mean_' if plot_type == 'gene_categories' else ''}gene_expr_"
                                                         f"{analysis_label}_{feature_space}_space")
            
def plot_cat_gp_score_gene_expr_summary(adata,
                                        cats,
                                        gps,
                                        plot_type,
                                        plot_category,
                                        plot_group,
                                        feature_space,
                                        suptitle,
                                        cat_title,
                                        save_fig,
                                        fig_name):
        
    if plot_category == cell_type_key:
        palette = seqfish_mouse_organogenesis_cell_type_colors
    else:
        palette = latent_cluster_colors 
    # Plot selected gene program latent scores
    if plot_type == "gene_categories":
        ncols = 6
        fig_width = 36
        wspace = 0.155
    elif plot_type == "individual_genes":
        ncols = 2 + adata.uns["n_top_genes"]
        fig_width = 12 + (6 * adata.uns["n_top_genes"])
        wspace = 0.3
    fig, axs = plt.subplots(nrows=len(gps), ncols=ncols, figsize=(fig_width, 6*len(gps)))
    if axs.ndim == 1:
        axs = axs.reshape(1, -1)

    title = fig.suptitle(t=suptitle,
                         x=0.55,
                         y=(1.1 if len(gps) == 1 else 0.94),
                         fontsize=20)
    for i, gp in enumerate(gps):
        if feature_space == "latent":
            sc.pl.umap(adata,
                       color=plot_category,
                       palette=palette,
                       groups=plot_group,
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
            sc.pl.spatial(adata=adata[adata.obs["sample"] == feature_space.split("_")[1]],
                          color=plot_category,
                          palette=palette,
                          groups=plot_group,
                          ax=axs[i, 0],
                          spot_size=0.03,
                          title=cat_title,
                          legend_loc="on data",
                          na_in_legend=False,
                          show=False)
            sc.pl.spatial(adata=adata[adata.obs["sample"] == feature_space.split("_")[1]],
                          color=gps[i],
                          color_map="RdBu",
                          spot_size=0.03,
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
                                      spot_size=0.03,
                                      title=f"Weighted mean gene expression \n {gene_category.replace('_', ' ')} ({adata.uns[f'{gp}_gene_category_importances'][j]:.2f})",
                                      show=False)                        
                    axs[i, j+2].xaxis.label.set_visible(False)
                    axs[i, j+2].yaxis.label.set_visible(False)
                else:
                    axs[i, j+2].set_visible(False)
        elif plot_type == "individual_genes":
            for j in range(len(adata.uns[f"{gp}_top_genes"])):
                if feature_space == "latent":
                    sc.pl.umap(adata,
                               color=adata.uns[f"{gp}_top_genes"][j],
                               color_map=("Blues" if adata.uns[f"{gp}_top_gene_signs"][j] == "+" else "Reds"),
                               ax=axs[i, 2+j],
                               legend_loc="on data",
                               na_in_legend=False,
                               title=f"{adata.uns[f'{gp}_top_genes'][j]}: {adata.uns[f'{gp}_top_gene_importances'][j]:.2f} ({adata.uns[f'{gp}_top_gene_entities'][j][0]}; {adata.uns[f'{gp}_top_gene_signs'][j]})",
                               show=False)
                elif "physical" in feature_space:
                    sc.pl.spatial(adata=adata[adata.obs["sample"] == feature_space.split("_")[1]],
                                  color=adata.uns[f"{gp}_top_genes"][j],
                                  color_map=("Blues" if adata.uns[f"{gp}_top_gene_signs"][j] == "+" else "Reds"),
                                  legend_loc="on data",
                                  na_in_legend=False,
                                  ax=axs[i, 2+j],
                                  # groups=cats[i],
                                  spot_size=0.03,
                                  title=f"{adata.uns[f'{gp}_top_genes'][j]}: {adata.uns[f'{gp}_top_gene_importances'][j]:.2f} ({adata.uns[f'{gp}_top_gene_entities'][j][0]}; {adata.uns[f'{gp}_top_gene_signs'][j]})",
                                  show=False)
                axs[i, 2+j].xaxis.label.set_visible(False)
                axs[i, 2+j].yaxis.label.set_visible(False)
            for k in range(len(adata.uns[f"{gp}_top_genes"]), ncols - 2):
                axs[i, 2+k].set_visible(False)            

    # Save and display plot
    plt.subplots_adjust(wspace=wspace, hspace=0.275)
    if save_fig:
        fig.savefig(f"{figure_folder_path}/{fig_name}.png",
                    bbox_extra_artists=(title,),
                    bbox_inches="tight")
    plt.show()
    
    
def add_sub_cell_type(row, cell_type):  
    if row[cell_type_key] == cell_type:
        return row[cell_type_key] + " cluster " + row[cell_type_latent_cluster_key]
    return row[cell_type_key]


def add_cell_type_latent_cluster_emphasis(row, comparison_latent_clusters):  
    for comparison_latent_cluster in comparison_latent_clusters:
        if row[cell_type_latent_cluster_key] == comparison_latent_cluster:
            return "-1"
    return row[cell_type_latent_cluster_key]


def store_top_gps_summary(model,
                          top_gps,
                          file_name):
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
    top_gps_summary_df.to_csv(f"{model_artifacts_folder_path}/{file_name}")
    return top_gps_summary_df