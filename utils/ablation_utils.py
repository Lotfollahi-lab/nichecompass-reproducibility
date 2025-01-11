import gc
import math
import os
import time
from datetime import datetime

import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from plottable import ColumnDefinition, Table
from plottable.cmap import normed_cmap
from plottable.formatters import tickcross
from plottable.plots import bar

from nichecompass.benchmarking import compute_cas, compute_clisis, compute_gcs, compute_mlami, compute_benchmarking_metrics
from nichecompass.utils import create_new_color_dict

    
def compute_ablation_points(df,
                            group_col,
                            metric_cols,
                            metric_cols_weights,
                            sort_metric_col="total_score"):
    """Compute ablation points and return df with added column."""
    # Compute ablation points and ranks for sorting
    for metric_col in metric_cols:
        metric_col_points = (df.groupby([group_col])[metric_col].mean().rank(ascending=True)
                           .rename(f"{metric_col}_points"))
        df = df.merge(metric_col_points, on=[group_col])

    if sort_metric_col == "total_score":
        df["total_score"] = 0
        for metric_col_weight, metric_col in zip(metric_cols_weights,
                                                 metric_cols):
            df["total_score"] += metric_col_weight * df[f"{metric_col}_points"]
            df["total_score_rank"] = df["total_score"].rank(ascending=False)
        
    df.sort_values(by=[f"{sort_metric_col}_rank", group_col],
                   inplace=True,
                   ascending=True)
    
    return df


def compute_metrics(artifact_folder_path,
                    dataset,
                    task,
                    timestamps,
                    cell_type_key,
                    batch_key,
                    spatial_key,
                    latent_key,
                    n_jobs=1,
                    metrics=[
                        "gcs",
                        "cas",
                        "clisis",
                        "cnmi",
                        "nasw",
                        "basw",
                        "bgc",
                        "bilisi"],
                    file_name="temp.csv"):
    """
    Compute spatial conservation, biological conservation, niche identification, and batch correction metrics to
    evaluate an ablation task, and return a dataframe with the computed metrics.
    """
    # Initialize metrics dicts
    benchmark_dict_acc = {"dataset": [],
                          "timestamp": []}
    for metric in metrics:
        if batch_key is None:
            if metric in ["basw", "bgc", "bilisi"]:
                continue
            else:
                benchmark_dict_acc[metric] = [] 
        else:
            benchmark_dict_acc[metric] = []
    
    compute_knng_flag = True
    # For each model compute metrics and append results to metrics dict
    for i, timestamp in enumerate(timestamps):
        benchmark_dict_acc["dataset"].append(dataset)
        benchmark_dict_acc["timestamp"].append(timestamp)
        print(f"Loading {dataset} model with timestamp {timestamp}.")
        try:
            adata = sc.read_h5ad(f"{artifact_folder_path}/{dataset}/models/{task}/{timestamp}/{dataset}_{task}.h5ad")
        except:
            print(f"Could not load adata with path '{artifact_folder_path}/{dataset}/models/{task}/{timestamp}/{dataset}_{task}.h5ad' ...")
            for key in benchmark_dict_acc.keys():
                if key in metrics:
                    benchmark_dict_acc[key].append(0.0)
            continue
        
        if (batch_key is None) & (not compute_knng_flag):
            print(knng_dict)
            # Load spatial knn graph from first iteration to avoid recomputation
            if "spatial_15knng_connectivities" in knng_dict.keys():
                adata.obsp[f"nichecompass_spatial_15knng_connectivities"] = knng_dict["spatial_15knng_connectivities"]
                adata.obsp[f"nichecompass_spatial_15knng_distances"] = knng_dict["spatial_15knng_distances"]
                adata.uns[f"nichecompass_spatial_15knng_n_neighbors"] = 15
            if "spatial_50knng_connectivities" in knng_dict.keys():
                adata.obsp[f"nichecompass_spatial_50knng_connectivities"] = knng_dict["spatial_50knng_connectivities"]
                adata.obsp[f"nichecompass_spatial_50knng_distances"] = knng_dict["spatial_50knng_distances"]
                adata.uns[f"nichecompass_spatial_50knng_n_neighbors"] = 50
            if "spatial_90knng_connectivities" in knng_dict.keys():
                adata.obsp[f"nichecompass_spatial_90knng_connectivities"] = knng_dict["spatial_90knng_connectivities"]
                adata.obsp[f"nichecompass_spatial_90knng_distances"] = knng_dict["spatial_90knng_distances"]
                adata.uns[f"nichecompass_spatial_90knng_n_neighbors"] = 90
        
        benchmark_dict = compute_benchmarking_metrics(
                adata=adata,
                metrics=metrics,
                cell_type_key=cell_type_key,
                batch_key=batch_key,
                spatial_key=spatial_key,
                latent_key=latent_key,
                n_jobs=n_jobs,
                seed=0,
                mlflow_experiment_id=None)
        
        for key, value in benchmark_dict.items():
            benchmark_dict_acc[key].append(value)
        
        benchmark_df = pd.DataFrame(benchmark_dict_acc)
        
        # Store results in temp file
        benchmark_df.to_csv(f"{artifact_folder_path}/{task}/{file_name.replace('.csv', '')}_metrics_temp.csv")
        
        if compute_knng_flag:
            # Store spatial knn graph from first iteration to avoid recomputation
            knng_dict = {}
            if f"nichecompass_spatial_15knng_connectivities" in adata.obsp:
                knng_dict["spatial_15knng_connectivities"] = adata.obsp[f"nichecompass_spatial_15knng_connectivities"]
            if f"nichecompass_spatial_50knng_connectivities" in adata.obsp:
                knng_dict["spatial_50knng_connectivities"] = adata.obsp[f"nichecompass_spatial_50knng_connectivities"]
            if f"nichecompass_spatial_90knng_connectivities" in adata.obsp:
                knng_dict["spatial_90knng_connectivities"] = adata.obsp[f"nichecompass_spatial_90knng_connectivities"]
            if f"nichecompass_spatial_15knng_distances" in adata.obsp:
                knng_dict["spatial_15knng_distances"] = adata.obsp[f"nichecompass_spatial_15knng_distances"]
            if f"nichecompass_spatial_50knng_distances" in adata.obsp:
                knng_dict["spatial_50knng_distances"] = adata.obsp[f"nichecompass_spatial_50knng_distances"]
            if f"nichecompass_spatial_90knng_distances" in adata.obsp:
                knng_dict["spatial_90knng_distances"] = adata.obsp[f"nichecompass_spatial_90knng_distances"]
            compute_knng_flag = False

    return benchmark_df

    
def get_loss_weights(row):  
    return f"lambda_edge_recon_{row['lambda_edge_recon_']}_+_lambda_gene_expr_recon_{row['lambda_gene_expr_recon_']}"


def get_gp_mask(row):  
    return f"{'fc_1_layer_decoder' if row['add_fc_gps_instead_of_gp_dict_gps'] else 'custom_gp_mask_target_genes_ratio_' + str(row['nichenet_keep_target_genes_ratio'])}"


# alternative version
def get_loss_weights_combination_alternative(row):  
    return f"""{'No Gene Expression Prediction & No Edge Reconstruction' if (row['lambda_edge_recon_'] == 0.0) & (row['lambda_gene_expr_recon_'] == 0.0) else (
                'Only Gene Expression Prediction' if (row['lambda_edge_recon_'] == 0.0) & (row['lambda_gene_expr_recon_'] != 0.0) else (
                'Only Edge Reconstruction' if (row['lambda_edge_recon_'] != 0.0) & (row['lambda_gene_expr_recon_'] == 0.0) else (
                'Gene Expression Prediction Prioritized' if (row['lambda_edge_recon_']/row['lambda_gene_expr_recon_'] < (500000 / 300)) else (
                'Edge Reconstruction Prioritized' if (row['lambda_edge_recon_']/row['lambda_gene_expr_recon_'] > (500000 / 300)) else 'Balanced Gene Expression Prediction & Edge Reconstruction'))))}"""


def get_loss_weights_combination(row):  
    return f"""{'Neither' if (row['lambda_edge_recon_'] == 0.0) & (row['lambda_gene_expr_recon_'] == 0.0) else (
                'Only Spatial Gene Expr Recon' if (row['lambda_edge_recon_'] == 0.0) & (row['lambda_gene_expr_recon_'] != 0.0) else (
                'Only Edge Recon' if (row['lambda_edge_recon_'] != 0.0) & (row['lambda_gene_expr_recon_'] == 0.0) else 'Spatial Gene Expr Recon & Edge Recon'))}"""


def plot_ablation_points(df,
                         ablation_col,
                         ablation_col_width,
                         group_col,
                         metric_cols,
                         metric_col_width=1.5,
                         show=True,
                         save_dir=None,
                         save_name="ablation_results.svg"):
    """"""
    groups = df[group_col].unique().tolist()
    df = df.pivot(index=[ablation_col], columns=[group_col, "score_type"], values="score")
    df.reset_index(inplace=True)
    df.columns = ['_'.join(col).strip("_") for col in df.columns.values]
    if len(groups) > 1:
        for i, group in enumerate(groups):
            df[f"Total Points ({i})"] = df[[col for col in list(df.columns) if group in col]].sum(axis=1)
        df["Total Points (All)"] = df[[col for col in list(df.columns) if "Total Points" not in col]].sum(axis=1)
        df.sort_values(by=["Total Points (All)"], inplace=True, ascending=False)
    else:
        df["Total Points"] = df[[col for col in list(df.columns)]].sum(axis=1)
        df.sort_values(by=["Total Points"], inplace=True, ascending=False)
    
    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)

    column_definitions = [
        ColumnDefinition(name=ablation_col,
                         title=ablation_col.replace("_", " ").title(),
                         width=ablation_col_width,
                         textprops={"ha": "left", "weight": "bold"})]

    aggregate_cols = [col for col in list(df.columns) if "Total" in col]

    for i, group in enumerate(groups):
        if len(groups) > 1:
            group_number_string = f"({i})"
        else:
            group_number_string = ""
        # Circles for the metric columns
        column_definitions += [
            ColumnDefinition(
                name=f"{group}_{col}",
                title=col.upper().replace("_", " ").replace("POINTS", "Points"),
                width=metric_col_width,
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25}},
                cmap=cmap_fn(df[f"{group}_{col}"]),
                group=f"Metric Points \n {group.replace('_', ' ').title()} {group_number_string}",
                border="left" if j == 0 else None,
                formatter="{:.2f}")
            for j, col in enumerate(metric_cols)]

        # Circles for the aggregate columns
        column_definitions += [
            ColumnDefinition(
                name=col,
                title=col,
                width=metric_col_width,
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25}},
                cmap=cmap_fn(df[col]),
                group="Aggregates",
                border="left" if j == 0 else None,
                formatter="{:.2f}")
            for j, col in enumerate(aggregate_cols)]
        
    # Allow to manipulate text post-hoc (in illustrator)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, ax = plt.subplots(figsize=(len(df.columns) * metric_col_width + ablation_col_width, 3 + 0.3 * len(df.columns)))
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
            index_col=ablation_col,
        ).autoset_fontcolors(colnames=df.columns)
    if show:
        plt.show()
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)        
        fig.savefig(os.path.join(save_dir, save_name), facecolor=ax.get_facecolor(), dpi=300)
    return tab


def plot_metrics_table(df,
                       ablation_col,
                       ablation_col_width,
                       group_col,
                       metric_cols,
                       metric_col_weights,
                       metric_col_titles=None,
                       metric_col_width=1.5,
                       plot_width=15,
                       plot_height=10,
                       group_label_dict={"starmap_plus_mouse_cns": "STARmap PLUS Mouse CNS",
                                         "xenium_human_breast_cancer": "Xenium Human Breast Cancer",
                                         "vizgen_merfish_human_ovarian_cancer": "Vizgen MERFISH Human Ovarian Cancer"},
                       show=True,
                       save_dir=None,
                       save_name="ablation_results.svg"):
    """"""
    if metric_col_titles is None:
        metric_col_titles = [col.upper().replace("_", " ").replace(" ", "\n") for col in metric_cols]
    
    groups = df[group_col].unique().tolist()
    df = df.pivot(index=[ablation_col], columns=[group_col, "score_type"], values="score")
    df.reset_index(inplace=True)
    df.columns = ['_'.join(col).strip("_") for col in df.columns.values]
    if len(groups) > 1:
        sorted_metrics_col_list = []
        for i, group in enumerate(groups): 
            sorted_group_metrics_col_list = sorted([col for col in list(df.columns) if group in col],
                                           key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
            df[f"Overall Score ({i})"] = np.average(df[sorted_group_metrics_col_list],
                                                    weights=metric_col_weights,
                                                    axis=1)
            sorted_metrics_col_list.extend(sorted_group_metrics_col_list)
        df = df[[ablation_col] + sorted_metrics_col_list + [f"Overall Score ({i})" for i in range(len(groups))]]
        df["Overall Score (All)"] = np.average(df[sorted_metrics_col_list],
                                               weights=metric_col_weights * len(groups),
                                               axis=1)
        df.sort_values(by=["Overall Score (All)"], inplace=True, ascending=False)
        
    else:
        sorted_metrics_col_list = sorted([col for col in list(df.columns) if any(col.endswith(metric) for metric in metric_cols)],
                                  key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
        df = df[[ablation_col] + sorted_metrics_col_list]
        df["Overall Score"] = np.average(df[sorted_metrics_col_list],
                                         weights=metric_col_weights,
                                         axis=1)
        df.sort_values(by=["Overall Score"], inplace=True, ascending=False)
    
    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)

    column_definitions = [
        ColumnDefinition(name=ablation_col,
                         title=ablation_col.replace("_", " ").title(),
                         width=ablation_col_width,
                         textprops={"ha": "left", "weight": "bold"})]

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
                group=f"Dataset Metrics \n {group_label_dict[group]} {group_number_string}",
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
            index_col=ablation_col,
        ).autoset_fontcolors(colnames=df.columns)
    if show:
        plt.show()
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)        
        fig.savefig(os.path.join(save_dir, save_name), facecolor=ax.get_facecolor(), dpi=300)
    return tab


def plot_metrics_boxplot(fig_title,
                         df,
                         group_col,
                         metric_cols,
                         figure_folder_path,
                         file_name,
                         metric_col_titles=None,
                         order=None,
                         plot_ratio_active_gps=False,
                         save_fig=False):
    """"""
    if metric_col_titles is None:
        metric_col_titles = metric_cols
    if plot_ratio_active_gps:
        fig, axes = plt.subplots(len(metric_cols) + 1, 1, sharey=True, figsize=(10, 20))
        sns.boxplot(data=df, ax=axes[2], x="ratio_active_gps", y=group_col)
        axes[2].set_title("Ratio of Active Gene Programs")
    else:
        fig, axes = plt.subplots(len(metric_cols), 1, sharey=True, figsize=(10, len(metric_cols) * 5))
    fig.suptitle(fig_title, fontsize=25)
    for i, (metric_col, metric_col_title) in enumerate(zip(metric_cols, metric_col_titles)):
        sns.boxplot(data=df, ax=axes[i], x=metric_col, y=group_col, order=order)
        axes[i].set_title(metric_col_title)
        axes[i].set_xlabel("score")
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.945,
                        wspace=0.25,
                        hspace=0.25)
    if save_fig:
        plt.savefig(f"{figure_folder_path}/{file_name}",
                    bbox_inches="tight")
        

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def visualize_latent_embeddings(artifact_folder_path,
                                plot_label,
                                dataset,
                                task,
                                timestamps,
                                cat_key,
                                sample_key,
                                groups=None,
                                spot_size=50.,
                                save_fig=True):
    """"""
    # Load arbitrary model to plot physical space
    adata = sc.read_h5ad(f"{artifact_folder_path}/{dataset}/models/{task}/{timestamps[0]}/{dataset}_{task}.h5ad")
    cat_colors = create_new_color_dict(
        adata=adata,
        cat_key=cat_key) 
    samples = adata.obs[sample_key].unique().tolist()
    
    ncols_samples = min(3, len(samples))
    ncols_models = min(3, len(timestamps))
    nrows = (2 + int(len(samples) / 3) +
             int(len(timestamps) / 3))
    # Create plot of cell type annotations in physical and latent spaces
    fig = plt.figure(figsize=(20, len(samples) + 20))
    title = fig.suptitle(t=f"{plot_label} in NicheCompass " \
                           "Latent and Physical Space",
                         y=0.96,
                         x=0.55,
                         fontsize=20)
    spec1 = gridspec.GridSpec(ncols=ncols_samples,
                              nrows=nrows,
                              width_ratios=[1] * ncols_samples,
                              height_ratios=[2] * nrows)
    spec2 = gridspec.GridSpec(ncols=ncols_models,
                              nrows=nrows,
                              width_ratios=[1] * ncols_models,
                              height_ratios=[2] * nrows)
    axs = []
    
    for i, sample in enumerate(samples):
        axs.append(fig.add_subplot(spec1[i]))
        sc.pl.spatial(adata=adata[adata.obs[sample_key] == sample],
                      color=[cat_key],
                      groups=groups,                  
                      palette=cat_colors,
                      spot_size=spot_size,
                      title=f"{plot_label} in \n Physical Space \n"
                            f"(Sample: {sample})",
                      ax=axs[i],
                      show=False)

    # Load all models and plot latent embeddings
    for j, timestamp in enumerate(timestamps):
        axs.append(fig.add_subplot(spec2[(math.ceil(len(samples)/ncols_samples) * ncols_models) + j]))
        adata = sc.read_h5ad(f"{artifact_folder_path}/{dataset}/models/{task}/{timestamp}/{dataset}_{task}.h5ad")
        sc.pl.umap(adata=adata,
                   color=[cat_key],
                   groups=groups,
                   palette=cat_colors,
                   title=f"{plot_label} in NicheCompass Latent Space",
                   legend_loc=None,
                   ax=axs[len(samples) + j],
                   show=False)
    
    # Create and position shared legend
    handles, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(handles,
                     labels,
                     loc="center left",
                     bbox_to_anchor=(0.98, 0.5))
    axs[0].get_legend().remove()

    # Adjust, save and display plot
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    if save_fig:
        fig.savefig(file_path,
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()
    
    
def visualize_niches(artifact_folder_path,
                     dataset,
                     task,
                     timestamps,
                     sample_key,
                     latent_key,
                     latent_leiden_resolution=0.2,
                     latent_cluster_key="nichecompass_latent_clusters",
                     spot_size=0.03,
                     groups=None,
                     save_fig=False,
                     save_folder_path="",
                     file_name="",
                     file_format="png"):
    for timestamp in timestamps:
        adata = sc.read_h5ad(f"{artifact_folder_path}/{dataset}/models/{task}/{timestamp}/{dataset}_{task}.h5ad")
        samples = adata.obs[sample_key].unique().tolist()

        # Compute latent Leiden clustering
        sc.tl.leiden(adata=adata,
                     resolution=latent_leiden_resolution,
                     key_added=latent_cluster_key,
                     neighbors_key=latent_key)
        
        latent_cluster_colors = create_new_color_dict(
            adata=adata,
            cat_key=latent_cluster_key)
        
        print(latent_cluster_colors)
        print(adata.obs[latent_cluster_key])

        fig = plt.figure(figsize=(12, 14))
        title = fig.suptitle(t=f"NicheCompass Latent Clusters " \
                               "in Latent and Physical Space",
                             y=0.96,
                             x=0.55,
                             fontsize=20)
        spec1 = gridspec.GridSpec(ncols=1,
                                  nrows=2,
                                  width_ratios=[1],
                                  height_ratios=[3, 2])
        spec2 = gridspec.GridSpec(ncols=len(samples),
                                  nrows=2,
                                  width_ratios=[1] * len(samples),
                                  height_ratios=[3, 2])
        axs = []
        axs.append(fig.add_subplot(spec1[0]))
        sc.pl.umap(adata=adata,
                   color=[latent_cluster_key],
                   groups=groups,
                   palette=latent_cluster_colors,
                   title=f"Latent Clusters in Latent Space",
                   ax=axs[0],
                   show=False)
        for idx, sample in enumerate(samples):
            axs.append(fig.add_subplot(spec2[len(samples) + idx]))
            sc.pl.spatial(adata=adata[adata.obs[sample_key] == sample],
                          color=[latent_cluster_key],
                          groups=groups,
                          palette=latent_cluster_colors,
                          spot_size=spot_size,
                          title=f"Latent Clusters in Physical Space \n"
                                f"(Sample: {sample})",
                          legend_loc=None,
                          ax=axs[idx+1],
                          show=False)

        # Create and position shared legend
        handles, labels = axs[0].get_legend_handles_labels()
        lgd = fig.legend(handles,
                         labels,
                         loc="center left",
                         bbox_to_anchor=(0.98, 0.5))
        axs[0].get_legend().remove()

        # Adjust, save and display plot
        plt.subplots_adjust(wspace=0.2, hspace=0.25)
        if save_fig:
            fig.savefig(f"{save_folder_path}/{file_name}_{timestamp}.{file_format}",
                        bbox_extra_artists=(lgd, title),
                        bbox_inches="tight")
        plt.show()
        gc.collect()
        

def scale_metric(metric,
                 reverse=True):
    scaled_metric = math.atan(metric) / (math.pi / 2)
    if reverse:
        scaled_metric = 1 - scaled_metric
    return scaled_metric