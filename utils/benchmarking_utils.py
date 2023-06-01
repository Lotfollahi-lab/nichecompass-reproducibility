import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scib

from nichecompass.benchmarking import compute_clisis, compute_cas


def compute_latent_space_comparison(dataset,
                                    run_number,
                                    srt_data_results_folder_path,
                                    cell_type_colors,
                                    dataset_title_string,
                                    cell_type_key,
                                    condition_key,
                                    figure_folder_path,
                                    cell_type_groups=None,
                                    spot_size=0.03,
                                    included_models=["NicheCompass",
                                                     "GraphST",
                                                     "scVI"],
                                    save_fig=True):
    fig = plt.figure(constrained_layout=True, figsize=(25, 15))
    suptitle = plt.suptitle(f"Sample Integration: Latent Space Comparison ({dataset_title_string})",
                            fontsize=35,
                            y=1.085)
    subfigs = fig.subfigures(nrows=2, ncols=1, hspace=0.15, wspace=0.05)
    subfigs[0].suptitle("Batch Annotation", fontsize=25, y=1.075)
    subfigs[1].suptitle("Cell Type Annotation", fontsize=25, y=1.075)
    axs_0 = subfigs[0].subplots(nrows=1, ncols=3)
    axs_1 = subfigs[1].subplots(nrows=1, ncols=3)
    
    # Load model-specific data
    for i, model in enumerate(["NicheCompass",
                               "GraphST",
                               "scVI"]):
        if model in included_models:
            adata = sc.read_h5ad(f"{srt_data_results_folder_path}/{dataset}_{model.lower()}_sample_integration_method_benchmarking.h5ad")
            
            # Get UMAP features from specified run
            adata.obsm["X_umap"] = adata.obsm[f"{model.lower()}_latent_run{run_number}_X_umap"]
            
            # Plot UMAP with batch annotations
            sc.pl.umap(adata,
                       color=[condition_key],
                       size=240000/len(adata),
                       ax=axs_0[i],
                       show=False,
                       legend_loc="right margin" if i == 0 else None)
            if i == 0:
                handles, labels = axs_0[i].get_legend_handles_labels()
                lgd_0 = fig.legend(handles,
                                   labels,
                                   fontsize=15,
                                   loc="center left",
                                   bbox_to_anchor=(1.01, 0.77))
                axs_0[i].get_legend().remove()

            # Plot UMAP with cell type annotations
            sc.pl.umap(adata,
                       color=[cell_type_key],
                       palette=cell_type_colors,
                       size=240000/len(adata),
                       ax=axs_1[i],
                       show=False,
                       legend_loc="right margin" if i == 0 else None)
            if i == 0:
                handles, labels = axs_1[i].get_legend_handles_labels()
                lgd_1 = fig.legend(handles,
                                   labels,
                                   fontsize=15,
                                   loc="center left",
                                   bbox_to_anchor=(1.01, 0.23))
                axs_1[i].get_legend().remove()

            del(adata)
        else:
            axs_0[i].grid(False)
            axs_0[i].set_xticks([])
            axs_0[i].set_yticks([])
            axs_0[i].set_xlabel("UMAP1")
            axs_0[i].set_ylabel("UMAP2")
            axs_1[i].grid(False)
            axs_1[i].set_xticks([])
            axs_1[i].set_yticks([])
            axs_1[i].set_xlabel("UMAP1")
            axs_1[i].set_ylabel("UMAP2")
        axs_0[i].set_title(model, fontsize=20, pad=10)
        axs_1[i].set_title(model, fontsize=20, pad=10)
        
    if save_fig:
        fig.savefig(f"{figure_folder_path}/method_comparison_{dataset}_run{run_number}_latent"
                    f"{'_' + cell_type_groups.replace(' ', '_').lower() if cell_type_groups is not None else ''}.svg",
                    bbox_inches="tight",
                    bbox_extra_artists=(suptitle, lgd_0, lgd_1),
                    format="svg")
    plt.show()

    
def compute_batch_integration_metrics(dataset,
                                      condition_key,
                                      cell_type_key,
                                      srt_data_results_folder_path,
                                      metric_artifacts_folder_path,
                                      spatial_key="spatial",
                                      latent_key="latent",
                                      included_models=["NicheCompass",
                                                       "GraphST",
                                                       "scVI"]):
    metrics_dict = {"Dataset": [],
                    "Model": [],
                    "Run": [],
                    "CAS": [],
                    "CLISIS": [],
                    "ASW": [],
                    "ILISI": []}
    for model in included_models:
        # Load model-specific data
        adata = sc.read_h5ad(f"{srt_data_results_folder_path}/{dataset}_{model.lower()}_sample_integration_method_benchmarking.h5ad")
        
        metrics_dict["Dataset"].append(dataset)
        metrics_dict["Model"].append(model)
        
        # Compute spatial nearest neighbor graph once, then use precomputed one
        spatial_knng_key = f"{model.lower()}_spatial"

        # Compute metrics per run
        for run_number in range(1, 11):
            metrics_dict["Run"].append(run_number)
            
            # Use precomputed latent nearest neighbor graph
            latent_knng_run_key = f"{model.lower()}_{latent_key}_run{run_number}"
            latent_run_key = f"{model.lower()}_{latent_key}_run{run_number}"

            # Spatial conservation metrics
            metrics_dict["CAS"].append(compute_cas(
                adata=adata,
                cell_type_key=cell_type_key,
                condition_key=condition_key,
                spatial_knng_key=spatial_knng_key,
                latent_knng_key=latent_knng_run_key,
                spatial_key=spatial_key,
                latent_key=latent_run_key))
            metrics_dict["CLISIS"].append(compute_clisis(
                adata=adata,
                cell_type_key=cell_type_key,
                condition_key=condition_key,
                spatial_knng_key=spatial_knng_key,
                latent_knng_key=latent_knng_run_key,
                spatial_key=spatial_key,
                latent_key=latent_run_key))

            # Batch correction metrics
            metrics_dict["ASW"].append(scib.me.silhouette_batch(
                adata=adata,
                batch_key=condition_key,
                label_key=cell_type_key,
                embed=latent_run_key))
            metrics_dict["ILISI"].append(scib.me.ilisi_graph(
                adata=adata,
                batch_key=condition_key,
                type_="embed",
                use_rep=latent_run_key))

        metric_df = pd.DataFrame(metrics_dict)

        # Store metrics to disk
        metric_df.to_csv(f"{metric_artifacts_folder_path}/method_comparison_metrics.csv")
            
            
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_batch_integration_results(df,
                                   metric_cols = ["cas", "clisis", "asw", "ilisi"],
                                   show=True,
                                   save_dir=None,
                                   save_name="batch_integration_results.svg"):
    datasets = df["dataset"].unique().tolist()
    df = df.pivot(index=["model_name", "spatially_aware"], columns=["dataset", "score_type"], values="score")
    df.reset_index(inplace=True)
    df.columns = ['_'.join(col).strip("_") for col in df.columns.values]
    if len(datasets) > 1:
        for i, dataset in enumerate(datasets):
            df[f"Total ({i})"] = df[[col for col in list(df.columns) if dataset in col]].mean(axis=1)
        df["Total (All)"] = df[list(set(list(df.columns)) - set(["model_name", "spatially_aware"]))].mean(axis=1)
        df.sort_values(by=["Total (All)"], inplace=True, ascending=False)
    else:
        df["Total"] = df[list(set(list(df.columns)) - set(["model_name"]))].mean(axis=1)
        df.sort_values(by=["Total"], inplace=True, ascending=False)
    df.rename(columns={"model_name": "Model"}, inplace=True)
    
    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)

    column_definitions = [
        ColumnDefinition(name="Model",
                         title="Model",
                         width=1.5,
                         textprops={"ha": "left", "weight": "bold"}),
        ColumnDefinition(name="spatially_aware",
                         title="Spatially \n Aware",
                         width=1.,
                         formatter=tickcross)]

    aggregate_cols = [col for col in list(df.columns) if "Total" in col]

    for i, dataset in enumerate(datasets):
        if len(datasets) > 1:
            dataset_number_string = f"({i})"
        else:
            dataset_number_string = ""
        # Circles for the metric columns
        column_definitions += [
            ColumnDefinition(
                name=f"{dataset}_{col}",
                title=col.upper(),
                width=1,
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25}},
                cmap=cmap_fn(df[f"{dataset}_{col}"]),
                group=f"Metrics \n {dataset.replace('_', ' ').title()} {dataset_number_string}",
                border="left" if j == 0 else None,
                formatter="{:.2f}")
            for j, col in enumerate(metric_cols)]

        # Bars for the aggregate columns
        column_definitions += [
            ColumnDefinition(
                name=col,
                title=col,
                width=1.5,
                plot_fn=bar,
                plot_kw={
                    "cmap": truncate_colormap(matplotlib.cm.YlOrRd, 0, 0.8),
                    "plot_bg_bar": False,
                    "annotate": True,
                    "height": 0.9,
                    "formatter": "{:.2f}",
                },
                group="Aggregate Scores",
                border="left" if j == 0 else None)
            for j, col in enumerate(aggregate_cols)]
        
    print(column_definitions)
        
    # Allow to manipulate text post-hoc (in illustrator)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, ax = plt.subplots(figsize=(len(df.columns) * 1., 3 + 0.3 * len(df.columns)))
        tab = Table(
            df,
            cell_kw={
                "linewidth": 0,
                "edgecolor": "k"},
            column_definitions=column_definitions,
            ax=ax,
            row_dividers=True,
            footer_divider=True,
            textprops={"fontsize": 10, "ha": "center"},
            row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
            col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
            column_border_kw={"linewidth": 1, "linestyle": "-"},
            index_col="Model",
        ).autoset_fontcolors(colnames=df.columns)
    if show:
        plt.show()
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)        
        fig.savefig(os.path.join(save_dir, save_name), facecolor=ax.get_facecolor(), dpi=300)
    return tab