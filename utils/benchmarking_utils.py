import os
import pickle

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


default_color_dict = {
    "0": "#66C5CC",
    "1": "#F6CF71",
    "2": "#F89C74",
    "3": "#DCB0F2",
    "4": "#87C55F",
    "5": "#9EB9F3",
    "6": "#FE88B1",
    "7": "#C9DB74",
    "8": "#8BE0A4",
    "9": "#B497E7",
    "10": "#D3B484",
    "11": "#B3B3B3",
    "12": "#276A8C", # Royal Blue
    "13": "#DAB6C4", # Pink
    "14": "#C38D9E", # Mauve-Pink
    "15": "#9D88A2", # Mauve
    "16": "#FF4D4D", # Light Red
    "17": "#9B4DCA", # Lavender-Purple
    "18": "#FF9CDA", # Bright Pink
    "19": "#FF69B4", # Hot Pink
    "20": "#FF00FF", # Magenta
    "21": "#DA70D6", # Orchid
    "22": "#BA55D3", # Medium Orchid
    "23": "#8A2BE2", # Blue Violet
    "24": "#9370DB", # Medium Purple
    "25": "#7B68EE", # Medium Slate Blue
    "26": "#4169E1", # Royal Blue
    "27": "#FF8C8C", # Salmon Pink
    "28": "#FFAA80", # Light Coral
    "29": "#48D1CC", # Medium Turquoise
    "30": "#40E0D0", # Turquoise
    "31": "#00FF00", # Lime
    "32": "#7FFF00", # Chartreuse
    "33": "#ADFF2F", # Green Yellow
    "34": "#32CD32", # Lime Green
    "35": "#228B22", # Forest Green
    "36": "#FFD8B8", # Peach
    "37": "#008080", # Teal
    "38": "#20B2AA", # Light Sea Green
    "39": "#00FFFF", # Cyan
    "40": "#00BFFF", # Deep Sky Blue
    "41": "#4169E1", # Royal Blue
    "42": "#0000CD", # Medium Blue
    "43": "#00008B", # Dark Blue
    "44": "#8B008B", # Dark Magenta
    "45": "#FF1493", # Deep Pink
    "46": "#FF4500", # Orange Red
    "47": "#006400", # Dark Green
    "48": "#FF6347", # Tomato
    "49": "#FF7F50", # Coral
    "50": "#CD5C5C", # Indian Red
    "51": "#B22222", # Fire Brick
    "52": "#FFB83F",  # Light Orange
    "53": "#8B0000", # Dark Red
    "54": "#D2691E", # Chocolate
    "55": "#A0522D", # Sienna
    "56": "#800000", # Maroon
    "57": "#808080", # Gray
    "58": "#A9A9A9", # Dark Gray
    "59": "#C0C0C0", # Silver
    "60": "#9DD84A",
    "61": "#F5F5F5", # White Smoke
    "62": "#F17171", # Light Red
    "63": "#000000", # Black
    "64": "#FF8C42", # Tangerine
    "65": "#F9A11F", # Bright Orange-Yellow
    "66": "#FACC15", # Golden Yellow
    "67": "#E2E062", # Pale Lime
    "68": "#BADE92", # Soft Lime
    "69": "#70C1B3", # Greenish-Blue
    "70": "#41B3A3", # Turquoise
    "71": "#5EAAA8", # Gray-Green
    "72": "#72B01D", # Chartreuse
    "73": "#9CD08F", # Light Green
    "74": "#8EBA43", # Olive Green
    "75": "#FAC8C3", # Light Pink
    "76": "#E27D60", # Dark Salmon
    "77": "#C38D9E", # Mauve-Pink
    "78": "#937D64", # Light Brown
    "79": "#B1C1CC", # Light Blue-Gray
    "80": "#88A0A8", # Gray-Blue-Green
    "81": "#4E598C", # Dark Blue-Purple
    "82": "#4B4E6D", # Dark Gray-Blue
    "83": "#8E9AAF", # Light Blue-Grey
    "84": "#C0D6DF", # Pale Blue-Grey
    "85": "#97C1A9", # Blue-Green
    "86": "#4C6E5D", # Dark Green
    "87": "#95B9C7", # Pale Blue-Green
    "88": "#C1D5E0", # Pale Gray-Blue
    "89": "#ECDB54", # Bright Yellow
    "90": "#E89B3B", # Bright Orange
    "91": "#CE5A57", # Deep Red
    "92": "#C3525A", # Dark Red
    "93": "#B85D8E", # Berry
    "94": "#7D5295", # Deep Purple
    "-1" : "#E1D9D1",
    "None" : "#E1D9D1"
}


def create_new_color_dict(
        adata,
        cat_key,
        color_palette="default",
        overwrite_color_dict={"-1" : "#E1D9D1"},
        skip_default_colors=0):
    """
    Create a dictionary of color hexcodes for a specified category.

    Parameters
    ----------
    adata:
        AnnData object.
    cat_key:
        Key in ´adata.obs´ where the categories are stored for which color
        hexcodes will be created.
    color_palette:
        Type of color palette.
    overwrite_color_dict:
        Dictionary with overwrite values that will take precedence over the
        automatically created dictionary.
    skip_default_colors:
        Number of colors to skip from the default color dict.

    Returns
    ----------
    new_color_dict:
        The color dictionary with a hexcode for each category.
    """
    new_categories = adata.obs[cat_key].unique().tolist()
    if color_palette == "cell_type_30":
        # https://github.com/scverse/scanpy/blob/master/scanpy/plotting/palettes.py#L40
        new_color_dict = {key: value for key, value in zip(
            new_categories,
            ["#023fa5",
             "#7d87b9",
             "#bec1d4",
             "#d6bcc0",
             "#bb7784",
             "#8e063b",
             "#4a6fe3",
             "#8595e1",
             "#b5bbe3",
             "#e6afb9",
             "#e07b91",
             "#d33f6a",
             "#11c638",
             "#8dd593",
             "#c6dec7",
             "#ead3c6",
             "#f0b98d",
             "#ef9708",
             "#0fcfc0",
             "#9cded6",
             "#d5eae7",
             "#f3e1eb",
             "#f6c4e1",
             "#f79cd4",
             '#7f7f7f',
             "#c7c7c7",
             "#1CE6FF",
             "#336600"])}
    elif color_palette == "cell_type_20":
        # https://github.com/vega/vega/wiki/Scales#scale-range-literals (some adjusted)
        new_color_dict = {key: value for key, value in zip(
            new_categories,
            ['#1f77b4',
             '#ff7f0e',
             '#279e68',
             '#d62728',
             '#aa40fc',
             '#8c564b',
             '#e377c2',
             '#b5bd61',
             '#17becf',
             '#aec7e8',
             '#ffbb78',
             '#98df8a',
             '#ff9896',
             '#c5b0d5',
             '#c49c94',
             '#f7b6d2',
             '#dbdb8d',
             '#9edae5',
             '#ad494a',
             '#8c6d31'])}
    elif color_palette == "cell_type_10":
        # scanpy vega10
        new_color_dict = {key: value for key, value in zip(
            new_categories,
            ['#7f7f7f',
             '#ff7f0e',
             '#279e68',
             '#e377c2',
             '#17becf',
             '#8c564b',
             '#d62728',
             '#1f77b4',
             '#b5bd61',
             '#aa40fc'])}
    elif color_palette == "batch":
        # sns.color_palette("colorblind").as_hex()
        new_color_dict = {key: value for key, value in zip(
            new_categories,
            ['#0173b2', '#d55e00', '#ece133', '#ca9161', '#fbafe4',
             '#949494', '#de8f05', '#029e73', '#cc78bc', '#56b4e9',
             '#F0F8FF', '#FAEBD7', '#00FFFF', '#7FFFD4', '#F0FFFF',
             '#F5F5DC', '#FFE4C4', '#000000', '#FFEBCD', '#0000FF',
             '#8A2BE2', '#A52A2A', '#DEB887', '#5F9EA0', '#7FFF00',
             '#D2691E', '#FF7F50', '#6495ED', '#FFF8DC', '#DC143C'])}
    elif color_palette == "default":
        new_color_dict = {key: value for key, value in zip(new_categories, list(default_color_dict.values())[skip_default_colors:])}
    for key, val in overwrite_color_dict.items():
        new_color_dict[key] = val
    return new_color_dict
    

def plot_simple_metrics_table(df,
                              model_col,
                              model_col_width,
                              group_col,
                              metric_cols,
                              metric_col_weights,
                              metric_col_titles=None,
                              metric_col_width=1.5,
                              aggregate_col_width=1.5,
                              plot_width=15,
                              plot_height=10,
                              group_label_dict={"seqfish_mouse_organogenesis_embryo2": "seqFISH \n Mouse Organogenesis",
                                                "nanostring_cosmx_human_nsclc_batch5": "nanoString CosMx \n Human NSCLC",
                                                "nanostring_cosmx_human_nsclc": "nanoString CosMx \n Human NSCLC",
                                                "nanostring_cosmx_human_nsclc_subsample_10pct": "nanoString CosMx \n Human NSCLC (10% Subsample)",
                                                "vizgen_merfish_mouse_liver": "MERFISH \n Mouse Liver",
                                                "slideseqv2_mouse_hippocampus": "Slide-seqV2 \n Mouse Hippocampus",
                                                "slideseqv2_mouse_hippocampus_subsample_25pct": "Slide-seqV2 \n Mouse Hippocampus (25% Subsample)",
                                                "sim1_1105genes_10000locs_strongincrements": "Simulation Data",
                                                "stereoseq_mouse_embryo": "StereoSeq Mouse Embryo",
                                                "starmap_mouse_mpfc": "STARmap Mouse MPFC"},
                              overall_score_col="Overall Score",
                       show=True,       
                       save_dir=None,
                       save_name=f"benchmarking_results.png"):
    """"""
    if metric_col_titles is None:
        metric_col_titles = [col.upper().replace("_", " ").replace(" ", "\n") for col in metric_cols]
    
    groups = df[group_col].unique().tolist()
    df = df.pivot(index=[model_col], columns=[group_col, "score_type"], values="score")
    df.reset_index(inplace=True)
    df.columns = ['_'.join(col).strip("_") for col in df.columns.values]
    df.columns = [col if overall_score_col not in col else overall_score_col for col in df.columns.values]
    if len(groups) > 1:
        sorted_metrics_col_list = []
        for i, group in enumerate(groups):
            sorted_group_metrics_col_list = sorted([col for col in list(df.columns) if ((group  in col) & (("subsample" in group) == ("subsample" in col)))],
                                           key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
            print(sorted_group_metrics_col_list)
            df[f"Overall Score ({i})"] = np.average(df[sorted_group_metrics_col_list],
                                                    weights=metric_col_weights,
                                                    axis=1)
            sorted_metrics_col_list.extend(sorted_group_metrics_col_list)
        df = df[[model_col] + sorted_metrics_col_list + [f"Overall Score ({i})" for i in range(len(groups))]]

        overall_sore_metrics_df = df[sorted_metrics_col_list]
        
        print(sorted_metrics_col_list)
        df["Overall Score (All)"] = np.average(overall_sore_metrics_df,
                                               weights=metric_col_weights * len(groups),
                                               axis=1)
        df.sort_values(by=["Overall Score (All)"], inplace=True, ascending=False)
        
    else:
        sorted_metrics_col_list = sorted([col for col in list(df.columns) if any(col.endswith(metric) for metric in metric_cols)],
                                  key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
        df = df[[model_col] + sorted_metrics_col_list + [overall_score_col]]
        df.sort_values(by=[overall_score_col], inplace=True, ascending=False)
    
    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)

    column_definitions = [
        ColumnDefinition(name=model_col,
                         title=model_col.replace("_", " ").title(),
                         width=model_col_width,
                         textprops={"ha": "left", "weight": "bold"}),
        #ColumnDefinition(name="spatially_aware",
        #         title="Spatially \n Aware",
        #         width=1.,
        #         formatter=tickcross)
    ]
    
    aggregate_cols = [col for col in list(df.columns) if overall_score_col in col]

    for i, group in enumerate(groups):
        if len(groups) > 1:
            group_number_string = f"({i})"
        else:
            group_number_string = "\nRaw Metrics"
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

        # Bars for the aggregate columns
        column_definitions += [
            ColumnDefinition(
                name=col,
                title=col.replace(" ", "\n"),
                width=aggregate_col_width,
                plot_fn=bar,
                plot_kw={
                    "cmap": cmap_fn(df[col]),
                    "plot_bg_bar": False,
                    "annotate": True,
                    "xlim": (0, 1),
                    "height": 0.9,
                    "formatter": "{:.3f}",
                    "textprops": {"fontsize": 12}
                },
                group="Scaled\nAggregates",
                border="left" if j == 0 else None)
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


def plot_metrics_table(df,
                       model_col,
                       model_col_width,
                       group_col,
                       metric_cols,
                       metric_col_weights,
                       metric_col_titles=None,
                       metric_col_width=1.5,
                       aggregate_col_width=1.5,
                       plot_width=15,
                       plot_height=10,
                       group_label_dict={"seqfish_mouse_organogenesis_embryo2": "seqFISH \n Mouse Organogenesis (100%)",
                                         "seqfish_mouse_organogenesis_subsample_50pct_embryo2": "seqFISH \n Mouse Organogenesis (50%)",
                                         "seqfish_mouse_organogenesis_subsample_25pct_embryo2": "seqFISH \n Mouse Organogenesis (25%)",
                                         "seqfish_mouse_organogenesis_subsample_10pct_embryo2": "seqFISH \n Mouse Organogenesis (10%)",
                                         "seqfish_mouse_organogenesis_subsample_5pct_embryo2": "seqFISH \n Mouse Organogenesis (5%)",
                                         "seqfish_mouse_organogenesis_subsample_1pct_embryo2": "seqFISH \n Mouse Organogenesis (1%)",
                                         "nanostring_cosmx_human_nsclc_batch5": "nanoString CosMx \n Human NSCLC (100%)",
                                         "nanostring_cosmx_human_nsclc_subsample_50pct_batch5": "nanoString CosMx \n Human NSCLC (50%)",
                                         "nanostring_cosmx_human_nsclc_subsample_25pct_batch5": "nanoString CosMx \n Human NSCLC (25%)",
                                         "nanostring_cosmx_human_nsclc_subsample_10pct_batch5": "nanoString CosMx \n Human NSCLC (10%)",
                                         "nanostring_cosmx_human_nsclc_subsample_5pct_batch5": "nanoString CosMx \n Human NSCLC (5%)",
                                         "nanostring_cosmx_human_nsclc_subsample_1pct_batch5": "nanoString CosMx \n Human NSCLC (1%)",
                                         "vizgen_merfish_mouse_liver": "MERFISH \n Mouse Liver (100%)",
                                         "vizgen_merfish_mouse_liver_subsample_50pct": "MERFISH \n Mouse Liver (50%)",
                                         "vizgen_merfish_mouse_liver_subsample_25pct": "MERFISH \n Mouse Liver (25%)",
                                         "vizgen_merfish_mouse_liver_subsample_10pct": "MERFISH \n Mouse Liver (10%)",
                                         "vizgen_merfish_mouse_liver_subsample_5pct": "MERFISH \n Mouse Liver (5%)",
                                         "vizgen_merfish_mouse_liver_subsample_1pct": "MERFISH \n Mouse Liver (1%)",
                                         "slideseqv2_mouse_hippocampus": "SlideSeqV2 \n Mouse Hippocampus (100%)",
                                         "slideseqv2_mouse_hippocampus_subsample_50pct": "SlideSeqV2 \n Mouse Hippocampus (50%)",
                                         "slideseqv2_mouse_hippocampus_subsample_25pct": "SlideSeqV2 \n Mouse Hippocampus (25%)",
                                         "slideseqv2_mouse_hippocampus_subsample_10pct": "SlideSeqV2 \n Mouse Hippocampus (10%)",
                                         "slideseqv2_mouse_hippocampus_subsample_5pct": "SlideSeqV2 \n Mouse Hippocampus (5%)",
                                         "slideseqv2_mouse_hippocampus_subsample_1pct": "SlideSeqV2 \n Mouse Hippocampus (1%)",
                                         "seqfish_mouse_organogenesis": "seqFISH \n Mouse Organogenesis (100%)",
                                         "seqfish_mouse_organogenesis_subsample_50pct": "seqFISH \n Mouse Organogenesis (50%)",
                                         "seqfish_mouse_organogenesis_subsample_25pct": "seqFISH \n Mouse Organogenesis (25%)",
                                         "seqfish_mouse_organogenesis_subsample_10pct": "seqFISH \n Mouse Organogenesis (10%)",
                                         "seqfish_mouse_organogenesis_subsample_5pct": "seqFISH \n Mouse Organogenesis (5%)",
                                         "seqfish_mouse_organogenesis_subsample_1pct": "seqFISH \n Mouse Organogenesis (1%)",
                                         "nanostring_cosmx_human_nsclc": "nanoString CosMx \n Human NSCLC (100%)",
                                         "nanostring_cosmx_human_nsclc_subsample_50pct": "nanoString CosMx \n Human NSCLC (50%)",
                                         "nanostring_cosmx_human_nsclc_subsample_25pct": "nanoString CosMx \n Human NSCLC (25%)",
                                         "nanostring_cosmx_human_nsclc_subsample_10pct": "nanoString CosMx \n Human NSCLC (10%)",
                                         "nanostring_cosmx_human_nsclc_subsample_5pct": "nanoString CosMx \n Human NSCLC (5%)",
                                         "nanostring_cosmx_human_nsclc_subsample_1pct": "nanoString CosMx \n Human NSCLC (1%)",
                                         "seqfish_mouse_organogenesis_imputed": "seqFISH \n Mouse Organogenesis (Imputed; 100%)",
                                         "seqfish_mouse_organogenesis_imputed_subsample_50pct": "seqFISH \n Mouse Organogenesis (Imputed; 50%)",
                                         "seqfish_mouse_organogenesis_imputed_subsample_25pct": "seqFISH \n Mouse Organogenesis (Imputed; 25%)",
                                         "seqfish_mouse_organogenesis_imputed_subsample_10pct": "seqFISH \n Mouse Organogenesis (Imputed; 10%)",
                                         "seqfish_mouse_organogenesis_imputed_subsample_5pct": "seqFISH \n Mouse Organogenesis (Imputed; 5%)",
                                         "seqfish_mouse_organogenesis_imputed_subsample_1pct": "seqFISH \n Mouse Organogenesis (Imputed; 1%)",},
                       show=True,
                       save_dir=None,
                       save_name=f"benchmarking_results.png"):
    """"""
    if metric_col_titles is None:
        metric_col_titles = [col.upper().replace("_", " ").replace(" ", "\n") for col in metric_cols]
    
    groups = df[group_col].unique().tolist()
    df = df.pivot(index=[model_col], columns=[group_col, "score_type"], values="score")
    df.reset_index(inplace=True)
    df.columns = ['_'.join(col).strip("_") for col in df.columns.values]
    if len(groups) > 1:
        sorted_metrics_col_list = []
        for i, group in enumerate(groups):
            sorted_group_metrics_col_list = sorted([col for col in list(df.columns) if (("Overall" not in col) & (group  in col) & (("subsample" in group) == ("subsample" in col)))],
                                           key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
            df.rename(columns={f"{group}_Overall Score": f"Overall Score ({i})"}, inplace=True)
            sorted_metrics_col_list.extend(sorted_group_metrics_col_list)
        df = df[[model_col] + sorted_metrics_col_list + [f"Overall Score ({i})" for i in range(len(groups))]]
        df.sort_values(by=[f"Overall Score (0)"], inplace=True, ascending=False)
        
        # Create separate dataframe to fill nan values with values of worst method
        overall_sore_metrics_df = df[sorted_metrics_col_list].copy()
        for column in overall_sore_metrics_df.columns:
            min_value = overall_sore_metrics_df[column].min()
            overall_sore_metrics_df.fillna({column: min_value}, inplace=True)
    else:
        sorted_metrics_col_list = sorted([col for col in list(df.columns) if any(col.endswith(metric) for metric in metric_cols)],
                                  key=lambda x: [metric_cols.index(metric) for metric in metric_cols if x.endswith(metric)])
        df = df[[model_col] + sorted_metrics_col_list]
        df["Overall Score"] = np.average(df[sorted_metrics_col_list],
                                         weights=metric_col_weights,
                                         axis=1)
        df.sort_values(by=["Overall Score"], inplace=True, ascending=False)
        df.drop("Overall Score (All)", axis=1, inplace=True)
    
    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)

    column_definitions = [
        ColumnDefinition(name=model_col,
                         title=model_col.replace("_", " ").title(),
                         width=model_col_width,
                         textprops={"ha": "left", "weight": "bold"}),
        #ColumnDefinition(name="spatially_aware",
        #         title="Spatially \n Aware",
        #         width=1.,
        #         formatter=tickcross)
    ]
    
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

        # Bars for the aggregate columns
        column_definitions += [
            ColumnDefinition(
                name=col,
                title=col.replace(" ", "\n"),
                width=aggregate_col_width,
                plot_fn=bar,
                plot_kw={
                    "cmap": cmap_fn(df[col]),
                    "plot_bg_bar": False,
                    "annotate": True,
                    "xlim": (0, 1),
                    "height": 0.9,
                    "formatter": "{:.3f}",
                    "textprops": {"fontsize": 12}
                },
                group="Aggregates",
                border="left" if j == 0 else None)
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


def plot_category_in_latent_and_physical_space(
        adata,
        plot_label,
        model_label,
        cat_key,
        groups,
        sample_key,
        samples,
        cat_colors,
        size,
        spot_size,
        save_fig,
        file_path,
        figsize=None):
    """Plot latent clusters in latent and physical space."""
    ncols = 1
    if samples is None:
        samples = []
    if figsize is None:
        figsize = (10, 20 + (len(samples) - 1) * 10)
    # Create plot of cell type annotations in physical and latent space
    fig = plt.figure(figsize=figsize)
    title = fig.suptitle(t=f"{plot_label} in {model_label} " \
                           "Latent and Physical Space",
                         y=0.96,
                         x=0.55,
                         fontsize=20)
    spec1 = gridspec.GridSpec(ncols=1,
                              nrows=1 + 1 + int(len(samples) / 3),
                              width_ratios=[1],
                              height_ratios=[3.5 - len(samples)/2] + [2] * (1 + int(len(samples) / 3)))
    spec2 = gridspec.GridSpec(ncols=ncols,
                              nrows=2 + len(samples),
                              width_ratios=[1] * ncols,
                              height_ratios=[2 + len(samples)] + [2] * (1 + len(samples)))
    axs = []
    axs.append(fig.add_subplot(spec1[0]))
    sc.pl.umap(adata=adata,
               color=[cat_key],
               groups=groups,
               palette=cat_colors,
               size=size,
               title=model_label,
               ax=axs[0],
               show=False)
    for idx, sample in enumerate(samples):
        axs.append(fig.add_subplot(spec2[ncols + idx]))
        sc.pl.spatial(adata=adata[adata.obs[sample_key] == sample],
                      color=[cat_key],
                      groups=groups,                  
                      palette=cat_colors,
                      spot_size=spot_size,
                      title="Physical Space",
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

    #for handle in lgd.legendHandles:
    #    handle.set_sizes([100])  # Adjust the size as needed
    axs[0].get_legend().remove()

    # Adjust, save and display plot
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    if save_fig:
        fig.savefig(file_path,
                    bbox_extra_artists=(lgd, title),
                    bbox_inches="tight")
    plt.show()
    

def compute_latent_space_comparison(dataset,
                                    run_number,
                                    benchmarking_folder_path,
                                    cell_type_colors,
                                    dataset_title_string,
                                    cell_type_key,
                                    condition_key,
                                    figure_folder_path,
                                    cell_type_groups=None,
                                    spot_size=0.03,
                                    included_models=["NicheCompass GCN",
                                                     "GraphST",
                                                     "scVI"],
                                    save_fig=True):
    fig = plt.figure(constrained_layout=True, figsize=(25, 15))
    suptitle = plt.suptitle(f"Latent Space Comparison ({dataset_title_string})",
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
            adata = sc.read_h5ad(f"{benchmarking_folder_path}/{dataset}_{model.lower()}.h5ad")
            
            # Get UMAP features from specified run
            adata.obsm["X_umap"] = adata.obsm[f"{model.lower()}_latent_run{run_number}_X_umap"]
            
            if condition_key is None:
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
        

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap