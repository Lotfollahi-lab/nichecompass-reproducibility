import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scib
import seaborn as sns

from nichecompass.benchmarking import compute_cas, compute_cca, compute_clisis, compute_gcs, compute_mlami


def compute_metrics(artifact_folder_path,
                    dataset,
                    task,
                    timestamps,
                    cell_type_key,
                    condition_key,
                    spatial_knng_key,
                    latent_knng_key,
                    spatial_key,
                    latent_key):
    """"""
    metrics_dict = {"dataset": [],
                    "timestamp": [],
                    "cca": [],
                    "cas": [],
                    #"clisis": []
                   }
    
    if condition_key is None:
        metrics_dict["gcs"] = []
        metrics_dict["mlami"] = []
    else:
        metrics_dict["asw"] = []
        metrics_dict["mlami"] = []        
    
    for timestamp in timestamps:
        # Load model-specific data
        adata = sc.read_h5ad(f"{artifact_folder_path}/{dataset}/models/{task}/{timestamp}/{dataset}_{task}.h5ad")
        
        metrics_dict["dataset"].append(dataset)
        metrics_dict["timestamp"].append(timestamp)
        
        # Cell identity conservation metrics
        metrics_dict["cca"].append(compute_cca(
                adata=adata,
                cell_cat_key=cell_type_key,
                latent_key=latent_key))
        
        # Spatial conservation metrics
        metrics_dict["cas"].append(compute_cas(
            adata=adata,
            cell_type_key=cell_type_key,
            condition_key=condition_key,
            spatial_knng_key=spatial_knng_key,
            latent_knng_key=latent_knng_key,
            spatial_key=spatial_key,
            latent_key=latent_key))
        """
        metrics_dict["clisis"].append(compute_clisis(
            adata=adata,
            cell_type_key=cell_type_key,
            condition_key=condition_key,
            spatial_knng_key=spatial_knng_key,
            latent_knng_key=latent_knng_key,
            spatial_key=spatial_key,
            latent_key=latent_key))
        """
        if condition_key is None:
            metrics_dict["gcs"].append(compute_gcs(
                adata=adata,
                spatial_knng_key=spatial_knng_key,
                latent_knng_key=latent_knng_key))   
            metrics_dict["mlami"].append(compute_mlami(
                adata=adata,
                spatial_knng_key=spatial_knng_key,
                latent_knng_key=latent_knng_key))
        else:
            # Batch correction metrics
            metrics_dict["asw"].append(scib.me.silhouette_batch(
                adata=adata,
                batch_key=condition_key,
                label_key=cell_type_key,
                embed=latent_key))
            metrics_dict["ilisi"].append(scib.me.ilisi_graph(
                adata=adata,
                batch_key=condition_key,
                type_="embed",
                use_rep=latent_key))
    metrics_df = pd.DataFrame(metrics_dict)
    return metrics_df
    
    
def get_loss_weights(row):  
    return f"lambda_edge_recon_{row['lambda_edge_recon_']}_+_lambda_gene_expr_recon_{row['lambda_gene_expr_recon_']}"
        

def plot_metrics(fig_title,
                 df,
                 group_col,
                 metric_cols,
                 metric_cols_weights,
                 sort_metric_col="total_score",
                 plot_ratio_active_gps=False,
                 save_fig=False,
                 file_name="ablation_metrics.png"):
    """"""
    # Compute evaluation metric ranks for sorting
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
    
    print(df)
    
    if plot_ratio_active_gps:
        fig, axes = plt.subplots(len(metric_cols) + 1, 1, sharey=True, figsize=(10, 20))
        sns.boxplot(data=df, ax=axes[2], x="ratio_active_gps", y=group_col)
        axes[2].set_title("Ratio of Active Gene Programs")
    else:
        fig, axes = plt.subplots(len(metric_cols), 1, sharey=True, figsize=(10, 15))
    fig.suptitle(fig_title, fontsize=15)
    for i, metric_col in enumerate(metric_cols):
        sns.boxplot(data=df, ax=axes[i], x=metric_col, y=group_col)
        axes[i].set_title(metric_col)
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.94,
                        wspace=0.175,
                        hspace=0.5)
    if save_fig:
        # Get time for timestamping saved artefacts
        now = datetime.now()
        current_timestamp = now.strftime("%d%m%Y_%H%M%S")
        benchmark_fig_run_dir = f"{figure_folder_path}/{dataset}/benchmarking/{experiment_name}/results/{current_timestamp}"
        os.makedirs(benchmark_fig_run_dir, exist_ok=True)
        plt.savefig(f"{benchmark_fig_run_dir}/{file_name}",
                    bbox_inches='tight')
        
    return df