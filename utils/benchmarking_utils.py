import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from plottable import ColumnDefinition, Table
from plottable.cmap import normed_cmap
from plottable.formatters import tickcross
from plottable.plots import bar

#from nichecompass.benchmarking import compute_benchmarking_metrics, compute_clisis, compute_cas


def plot_metrics_table(df,
                       model_col,
                       model_col_width,
                       group_col,
                       metric_cols,
                       metric_col_weights,
                       metric_col_titles=None,
                       metric_col_width=1.5,
                       plot_width=15,
                       plot_height=10,
                       group_label_dict={"seqfish_mouse_organogenesis_embryo2": "seqFISH \n Mouse Organogenesis \n (100%)",
                                         "seqfish_mouse_organogenesis_subsample_50pct_embryo2": "seqFISH \n Mouse Organogenesis (50%)",
                                         "seqfish_mouse_organogenesis_subsample_25pct_embryo2": "seqFISH \n Mouse Organogenesis (25%)",
                                         "seqfish_mouse_organogenesis_subsample_10pct_embryo2": "seqFISH \n Mouse Organogenesis (10%)",
                                         "seqfish_mouse_organogenesis_subsample_5pct_embryo2": "seqFISH \n Mouse Organogenesis (5%)",
                                         "seqfish_mouse_organogenesis_subsample_1pct_embryo2": "seqFISH \n Mouse Organogenesis (1%)",
                                         "nanostring_cosmx_human_nsclc_batch5": "nanoString CosMx \n Human NSCLC \n (100%)",
                                         "nanostring_cosmx_human_nsclc_subsample_50pct_batch5": "nanoString CosMx \n Human NSCLC (50%)",
                                         "nanostring_cosmx_human_nsclc_subsample_25pct_batch5": "nanoString CosMx \n Human NSCLC (25%)",
                                         "nanostring_cosmx_human_nsclc_subsample_10pct_batch5": "nanoString CosMx \n Human NSCLC (10%)",
                                         "nanostring_cosmx_human_nsclc_subsample_5pct_batch5": "nanoString CosMx \n Human NSCLC (5%)",
                                         "nanostring_cosmx_human_nsclc_subsample_1pct_batch5": "nanoString CosMx \n Human NSCLC (1%)",
                                         "vizgen_merfish_mouse_liver": "MERFISH \n Mouse Liver \n (100%)",
                                         "vizgen_merfish_mouse_liver_subsample_50pct": "MERFISH \n Mouse Liver (50%)",
                                         "vizgen_merfish_mouse_liver_subsample_25pct": "MERFISH \n Mouse Liver (25%)",
                                         "vizgen_merfish_mouse_liver_subsample_10pct": "MERFISH \n Mouse Liver (10%)",
                                         "vizgen_merfish_mouse_liver_subsample_5pct": "MERFISH \n Mouse Liver (5%)",
                                         "vizgen_merfish_mouse_liver_subsample_1pct": "MERFISH \n Mouse Liver (1%)",
                                         "slideseqv2_mouse_hippocampus": "SlideSeqV2 \n Mouse Hippocampus \n (100%)",
                                         "slideseqv2_mouse_hippocampus_subsample_50pct": "SlideSeqV2 \n Mouse Hippocampus (50%)",
                                         "slideseqv2_mouse_hippocampus_subsample_25pct": "SlideSeqV2 \n Mouse Hippocampus (25%)",
                                         "slideseqv2_mouse_hippocampus_subsample_10pct": "SlideSeqV2 \n Mouse Hippocampus (10%)",
                                         "slideseqv2_mouse_hippocampus_subsample_5pct": "SlideSeqV2 \n Mouse Hippocampus (5%)",
                                         "slideseqv2_mouse_hippocampus_subsample_1pct": "SlideSeqV2 \n Mouse Hippocampus (1%)"},
                       show=True,
                       save_dir=None,
                       save_name=f"benchmarking_results.png"):
    """"""
    if metric_col_titles is None:
        metric_col_titles = [col.upper().replace("_", " ").replace(" ", "\n") for col in metric_cols]
    
    groups = df[group_col].unique().tolist()
    df = df.pivot(index=[model_col, "spatially_aware"], columns=[group_col, "score_type"], values="score")
    df.reset_index(inplace=True)
    df.columns = ['_'.join(col).strip("_") for col in df.columns.values]
    if len(groups) > 1:
        sorted_metrics_col_list = []
        for i, group in enumerate(groups):
            sorted_group_metrics_col_list = sorted([col for col in list(df.columns) if ((group  in col) & (("subsample" in group) == ("subsample" in col)))],
                                           key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
            df[f"Overall Score ({i})"] = np.average(df[sorted_group_metrics_col_list],
                                                    weights=metric_col_weights,
                                                    axis=1)
            sorted_metrics_col_list.extend(sorted_group_metrics_col_list)
        df = df[[model_col, "spatially_aware"] + sorted_metrics_col_list + [f"Overall Score ({i})" for i in range(len(groups))]]
        df["Overall Score (All)"] = np.average(df[sorted_metrics_col_list],
                                               weights=metric_col_weights * len(groups),
                                               axis=1)
        df.sort_values(by=["Overall Score (All)"], inplace=True, ascending=False)
        
    else:
        sorted_metrics_col_list = sorted([col for col in list(df.columns) if any(col.endswith(metric) for metric in metric_cols)],
                                  key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
        df = df[[model_col, "spatially_aware"] + sorted_metrics_col_list]
        df["Overall Score"] = np.average(df[sorted_metrics_col_list],
                                         weights=metric_col_weights,
                                         axis=1)
        df.sort_values(by=["Overall Score"], inplace=True, ascending=False)
    
    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)

    column_definitions = [
        ColumnDefinition(name=model_col,
                         title=model_col.replace("_", " ").title(),
                         width=model_col_width,
                         textprops={"ha": "left", "weight": "bold"}),
        ColumnDefinition(name="spatially_aware",
                 title="Spatially \n Aware",
                 width=1.,
                 formatter=tickcross)]
    
    aggregate_cols = [col for col in list(df.columns) if "Overall" in col]

    for i, group in enumerate(groups):
        if len(groups) > 1:
            group_number_string = f"({i})"
        else:
            group_number_string = ""
        # Circles for the metric columns
        column_definitions += [
            ColumnDefinition(
                name=f"{group}_{col}",
                title=metric_col_titles[j],
                width=metric_col_width,
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25}},
                cmap=cmap_fn(df[f"{group}_{col}"]),
                group=f"{group_label_dict[group]} {group_number_string}",
                border="left" if j == 0 else None,
                formatter="{:.3f}")
            for j, col in enumerate(metric_cols)]

        # Circles for the aggregate columns
        column_definitions += [
            ColumnDefinition(
                name=col,
                title=col.replace(" ", "\n"),
                width=metric_col_width,
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25}},
                cmap=cmap_fn(df[col]),
                group="Aggregates",
                border="left" if j == 0 else None,
                formatter="{:.3f}")
            for j, col in enumerate(aggregate_cols)]
        
    # Allow to manipulate text post-hoc (in illustrator)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, ax = plt.subplots(figsize=(plot_width, plot_height))
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
            index_col=model_col,
        ).autoset_fontcolors(colnames=df.columns)
    if show:
        plt.show()
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)        
        fig.savefig(os.path.join(save_dir, save_name), facecolor=ax.get_facecolor(), dpi=300)
    return tab


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


def compute_combined_benchmarking_metrics(model_adata,
                                          model_name,
                                          metrics,
                                          cell_type_key="cell_type",
                                          spatial_key="spatial",
                                          run_number_list=list(np.arange(1, 9)),
                                          n_neighbors_list=[4, 4, 8, 8, 12, 12, 16, 16]):
    benchmarking_dict_list = []
    for run_number, n_neighbors in zip(run_number_list, n_neighbors_list):
        
        print(model_adata)
        
        # Compute NicheCompass benchmarking metrics
        benchmarking_dict = compute_benchmarking_metrics(adata=model_adata,
                                                         metrics=metrics,
                                                         cell_type_key=cell_type_key,
                                                         spatial_key=spatial_key,
                                                         latent_key=f"{model_name}_latent_run{run_number}")
        benchmarking_dict["model_name"] = model_name
        benchmarking_dict["run"] = run_number
        benchmarking_dict_list.append(benchmarking_dict)
    return benchmarking_dict_list


def compute_combined_benchmarking_metrics_for_all_models(
        dataset,
        metrics=["gcs",
                 "mlami",
                 "cas",
                 "clisis",
                 "cnmi",
                 "casw",
                 "clisi",
                 "nasw",
                 "basw",
                 "bgc",
                 "blisi",
                 "kbet"],
        cell_type_key="cell_type",
        run_number_list=list(np.arange(1, 9)),
        n_neighbors_list=[4, 4, 8, 8, 12, 12, 16, 16],
        included_models=["nichecompass",
                         "deeplinc",
                         "graphst",
                         "sagenet",
                         "expimap",
                         "scvi"],
        benchmarking_type="single_sample_method_benchmarking"):
    # Configure dataset artifact folder path
    benchmarking_folder_path = f"../../artifacts/{benchmarking_type}/"
    os.makedirs(benchmarking_folder_path, exist_ok=True)
    
    benchmarking_dict_list = []
    
    # scVI
    if "scvi" in included_models:
        print("Computing metrics for scVI...")
        adata_scvi = sc.read_h5ad(benchmarking_folder_path + f"{dataset}_scvi.h5ad")
        benchmarking_dict_list_scvi = compute_combined_benchmarking_metrics(
            model_adata=adata_scvi,
            model_name="scvi",
            metrics=metrics,
            run_number_list=run_number_list,
            n_neighbors_list=n_neighbors_list,
            cell_type_key=cell_type_key)  
        benchmarking_dict_list += benchmarking_dict_list_scvi
        with open(f"{benchmarking_folder_path}/benchmarking_dict_list.pickle", "wb") as f:
            pickle.dump(benchmarking_dict_list, f)
        del(adata_scvi)
        print("")

    # expiMap
    if "expimap" in included_models:
        print("Computing metrics for expiMap...")
        adata_expimap = sc.read_h5ad(benchmarking_folder_path + f"{dataset}_expimap.h5ad")
        benchmarking_dict_list_expimap = compute_combined_benchmarking_metrics(
            model_adata=adata_expimap,
            model_name="expimap",
            metrics=metrics,
            run_number_list=run_number_list,
            n_neighbors_list=n_neighbors_list,
            cell_type_key=cell_type_key)  
        benchmarking_dict_list += benchmarking_dict_list_expimap
        with open(f"{benchmarking_folder_path}/benchmarking_dict_list.pickle", "wb") as f:
            pickle.dump(benchmarking_dict_list, f)
        del(adata_expimap)
        print("")
    
    # SageNet
    if "sagenet" in included_models:
        print("Computing metrics for SageNet...")
        adata_sagenet = sc.read_h5ad(benchmarking_folder_path + f"{dataset}_sagenet.h5ad")
        benchmarking_dict_list_sagenet = compute_combined_benchmarking_metrics(
            model_adata=adata_sagenet,
            model_name="sagenet",
            metrics=metrics,
            run_number_list=run_number_list,
            n_neighbors_list=n_neighbors_list,
            cell_type_key=cell_type_key) 
        benchmarking_dict_list += benchmarking_dict_list_sagenet
        with open(f"{benchmarking_folder_path}/benchmarking_dict_list.pickle", "wb") as f:
            pickle.dump(benchmarking_dict_list, f)
        del(adata_sagenet)
        print("")
    
    # DeepLinc
    if "deeplinc" in included_models:
        print("Computing metrics for DeepLinc...")
        adata_deeplinc = sc.read_h5ad(benchmarking_folder_path + f"{dataset}_deeplinc.h5ad")
        benchmarking_dict_list_deeplinc = compute_combined_benchmarking_metrics(
            model_adata=adata_deeplinc,
            model_name="deeplinc",
            metrics=metrics,
            run_number_list=run_number_list,
            n_neighbors_list=n_neighbors_list,
            cell_type_key=cell_type_key)
        benchmarking_dict_list += benchmarking_dict_list_deeplinc
        with open(f"{benchmarking_folder_path}/benchmarking_dict_list.pickle", "wb") as f:
            pickle.dump(benchmarking_dict_list, f)
        del(adata_deeplinc)
        print("")
    
    # GraphST
    if "graphst" in included_models:
        print("Computing metrics for GraphST...")
        adata_graphst = sc.read_h5ad(benchmarking_folder_path + f"{dataset}_graphst.h5ad")
        benchmarking_dict_list_graphst = compute_combined_benchmarking_metrics(
            model_adata=adata_graphst,
            model_name="graphst",
            metrics=metrics,
            run_number_list=run_number_list,
            n_neighbors_list=n_neighbors_list,
            cell_type_key=cell_type_key)
        benchmarking_dict_list += benchmarking_dict_list_graphst
        with open(f"{benchmarking_folder_path}/benchmarking_dict_list.pickle", "wb") as f:
            pickle.dump(benchmarking_dict_list, f)
        del(adata_graphst)
        print("")

    # NicheCompass
    if "nichecompass" in included_models:
        print("Computing metrics for NicheCompass...")
        adata_nichecompass = sc.read_h5ad(benchmarking_folder_path + f"{dataset}_nichecompass.h5ad")
        benchmarking_dict_list_nichecompass = compute_combined_benchmarking_metrics(
            model_adata=adata_nichecompass,
            model_name="nichecompass",
            metrics=metrics,
            run_number_list=run_number_list,
            n_neighbors_list=n_neighbors_list,
            cell_type_key=cell_type_key)
        benchmarking_dict_list += benchmarking_dict_list_nichecompass
        with open(f"{benchmarking_folder_path}/benchmarking_dict_list.pickle", "wb") as f:
            pickle.dump(benchmarking_dict_list, f)
        del(adata_nichecompass)
        print("")
        

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_benchmarking_results(df,
                              show=True,
                              save_dir=None,
                              save_name="benchmarking_results.svg"):
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

    metric_cols = ["cas", "clisis", "gcs", "mlami"]
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
                width=1,
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