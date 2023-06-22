default_color_dict = {
    "0": "#FFD8B8",  # Peach
    "1": "#FFAA80",  # Light Coral
    "2": "#FF4D4D",  # Light Red
    "3": "#FF8C8C",  # Salmon Pink
    "4": "#FFD966",  # Yellow
    "5": "#FFB83F",  # Light Orange
    "6": "#FF9326",  # Orange
    "7": "#FFED80",  # Soft Yellow
    "8": "#FFFEA3",  # Pale Canary
    "9": "#A5DEE4",  # Pale Blue
    "10": "#77B1BD", # Sky Blue
    "11": "#4E8CA7", # Light Blue
    "12": "#276A8C", # Royal Blue
    "13": "#DAB6C4", # Pink
    "14": "#C38D9E", # Mauve-Pink
    "15": "#9D88A2", # Mauve
    "16": "#6B5B95", # Purple
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
    "27": "#4682B4", # Steel Blue
    "28": "#00CED1", # Dark Turquoise
    "29": "#48D1CC", # Medium Turquoise
    "30": "#40E0D0", # Turquoise
    "31": "#00FF00", # Lime
    "32": "#7FFF00", # Chartreuse
    "33": "#ADFF2F", # Green Yellow
    "34": "#32CD32", # Lime Green
    "35": "#228B22", # Forest Green
    "36": "#006400", # Dark Green
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
    "47": "#FF0000", # Red
    "48": "#FF6347", # Tomato
    "49": "#FF7F50", # Coral
    "50": "#CD5C5C", # Indian Red
    "51": "#B22222", # Fire Brick
    "52": "#A52A2A", # Brown
    "53": "#8B0000", # Dark Red
    "54": "#D2691E", # Chocolate
    "55": "#A0522D", # Sienna
    "56": "#800000", # Maroon
    "57": "#808080", # Gray
    "58": "#A9A9A9", # Dark Gray
    "59": "#C0C0C0", # Silver
    "60": "#D3D3D3", # Light Gray
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
    "-1" : "#9DD84A",
    "None" : "#D3D3D3"
}

def create_new_color_dict(
        adata,
        cat_key):
    new_categories = adata.obs[cat_key].unique().tolist()
    new_color_dict = {key: value for key, value in zip(new_categories, default_color_dict.values())}
    return new_color_dict
    
latent_cluster_colors = {
    "0": "#FFD8B8", # Peach
    "1": "#F9CB9C", # Light Peach
    "2": "#FFB8A2", # Apricot
    "3": "#F08CAE", # Pink
    "4": "#F49FAD", # Pale Pink
    "5": "#E7F0A4", # Pale Canary
    "6": "#C7F2C2", # Light Green
    "7": "#97E6A1", # Soft Green
    "8": "#FE938C", # Coral
    "9": "#A5DEE4", # Pale Blue
    "10": "#5CA4A9", # Blue-Green
    "11": "#6B5B95", # Lavender
    "12": "#F9A03F", # Orange
    "13": "#F7DB6A", # Light Yellow
    "14": "#EEBB4D", # Light Amber
    "15": "#D6E4B2", # Pale Green
    "16": "#A8DADC", # Pale Cyan
    "17": "#3D5A80", # Dark Blue
    "18": "#3E3F8A", # Navy Blue
    "19": "#218380", # Teal
    "20": "#90BE6D", # Soft Green-Yellow
    "21": "#FFD369", # Yellow
    "22": "#ED553B", # Red-Orange
    "23": "#DA627D", # Mauve
    "24": "#6C5B7B", # Purple
    "25": "#4ECDC4", # Mint
    "26": "#65AADD", # Sky Blue
    "27": "#8FBFE0", # Powder Blue
    "28": "#A2D2FF", # Pale Sky Blue
    "29": "#F3C969", # Light Amber-Yellow
    "30": "#EE6C4D", # Light Red-Orange
    "31": "#EC4E20", # Bright Red
    "32": "#D64161", # Dark Pink
    "33": "#FF7A5A", # Bright Coral
    "34": "#E7A977", # Light Coral
    "35": "#FECE44", # Bright Yellow
    "36": "#FFC55F", # Yellow-Orange
    "37": "#F89E7B", # Light Coral-Orange
    "38": "#7ECEFD", # Baby Blue
    "39": "#C9B1BD", # Pale Mauve
    "40": "#E6A0C4", # Light Pink
    "41": "#E36BAE", # Bright Pink
    "42": "#8B5B6E", # Mauve-Brown
    "43": "#748CAB", # Blue-Gray
    "44": "#E5E5E5", # Light Gray
    "45": "#C4C4C4", # Gray
    "46": "#A4A4A4", # Dark Gray
    "47": "#4D4D4D", # Charcoal
    "48": "#F8B195", # Dusty Peach
    "49": "#F67280", # Salmon Pink
    "50": "#C06C84", # Rose
    "51": "#6C5B7B", # Muted Purple
    "52": "#355C7D", # Dark Slate Blue
    "53": "#6C7B95", # Gray-Blue
    "54": "#D6BCC0", # Light Mauve
    "55": "#D5B9B2", # Dusty Rose
    "56": "#A56C7B", # Muted Rose
    "57": "#F4A261", # Light Orange
    "58": "#F29E4C", # Bright Orange
    "59": "#E76F51", # Coral-Red
    "60": "#DA627D", # Mauve
    "61": "#9B4DCA", # Lavender-Purple
    "62": "#8A5B5C", # Brownish Red
    "63": "#9A8478", # Beige
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
    "-1" : "#9DD84A",
    "None" : "#D3D3D3"
}

starmap_plus_mouse_cns_cell_type_colors = {
    "Vascular and leptomeningeal cells": "#FFD8B8", # Peach
    "Unannotated": "#F9CB9C", # Light Peach
    "Pericytes": "#FFB8A2", # Apricot
    "Astrocytes": "#F08CAE", # Pink
    "Vascular smooth muscle cells": "#F49FAD", # Pale Pink
    "Oligodendrocytes": "#E7F0A4", # Pale Canary
    "Vascular endothelial cells": "#C7F2C2", # Light Green
    "Microglia": "#97E6A1", # Soft Green
    "Oligodendrocyte precursor cells": "#FE938C", # Coral
    "Olfactory ensheathing cells": "#A5DEE4", # Pale Blue
    "Telencephalon inhibitory interneurons": "#5CA4A9", # Blue-Green
    "Telencephalon projecting excitatory neurons": "#6B5B95", # Lavender
    "Non-glutamatergic neuroblasts": "#F9A03F", # Orange
    "Cholinergic and monoaminergic neurons": "#F7DB6A", # Light Yellow
    "Perivascular macrophages": "#EEBB4D", # Light Amber
    "Choroid plexus epithelial cells": "#D6E4B2", # Pale Green
    "Di- and mesencephalon inhibitory neurons": "#A8DADC", # Pale Cyan
    "Hindbrain neurons/Spinal cord neurons": "#3D5A80", # Dark Blue
    "Telencephalon projecting inhibitory neurons": "#3E3F8A", # Navy Blue
    "Olfactory inhibitory neurons": "#218380", # Teal
    "Di- and mesencephalon excitatory neurons": "#90BE6D", # Soft Green-Yellow
    "Glutamatergic neuroblasts": "#FFD369", # Yellow
    "Cerebellum neurons": "#ED553B", # Red-Orange
    "Peptidergic neurons": "#DA627D", # Mauve
    "Ependymal cells": "#6C5B7B", # Purple
    "Dentate gyrus granule neurons": "#4ECDC4", # Mint
    "Subcommissural organ hypendymal cells": "#65AADD", # Sky Blue 
}

seqfish_mouse_organogenesis_cell_type_colors = {
    "Epiblast" : "#FFD8B8", # Peach
    "Primitive Streak" : "#F9CB9C", # Light Peach
    "Caudal epiblast" : "#FFB8A2", # Apricot
    "PGC" : "#F08CAE", # Pink
    "Anterior Primitive Streak" : "#F49FAD", # Pale Pink
    "Notochord" : "#E7F0A4", # Pale Canary
    "Def. endoderm" : "#C7F2C2", # Light Green
    "Definitive endoderm" : "#97E6A1", # Soft Green
    "Gut" : "#FE938C", # Coral
    "Gut tube" : "#A5DEE4", # Pale Blue
    "Nascent mesoderm" : "#5CA4A9", # Blue-Green
    "Mixed mesoderm" : "#6B5B95", # Lavender
    "Intermediate mesoderm" : "#F9A03F", # Orange
    "Caudal Mesoderm" : "#F7DB6A", # Light Yellow
    "Paraxial mesoderm" : "#EEBB4D", # Light Amber
    "Somitic mesoderm" : "#D6E4B2", # Pale Green
    "Pharyngeal mesoderm" : "#A8DADC", # Pale Cyan
    "Splanchnic mesoderm" : "#3D5A80", # Dark Blue
    "Cardiomyocytes" : "#3E3F8A", # Navy Blue
    "Allantois" : "#218380", # Teal
    "ExE mesoderm" : "#90BE6D", # Soft Green-Yellow
    "Lateral plate mesoderm" : "#FFD369", # Yellow
    "Mesenchyme" : "#ED553B", # Red-Orange
    "Mixed mesenchymal mesoderm" : "#DA627D", # Mauve
    "Haematoendothelial progenitors" : "#6C5B7B", # Purple
    "Endothelium" : "#4ECDC4", # Mint
    "Blood progenitors 1" : "#65AADD", # Sky Blue
    "Blood progenitors 2" : "#8FBFE0", # Powder Blue
    "Erythroid1" : "#A2D2FF", # Pale Sky Blue
    "Erythroid2" : "#F3C969", # Light Amber-Yellow
    "Erythroid3" : "#EE6C4D", # Light Red-Orange
    "Erythroid" : "#EC4E20", # Bright Red
    "Blood progenitors" : "#D64161", # Dark Pink
    "NMP" : "#FF7A5A", # Bright Coral
    "Rostral neurectoderm" : "#E7A977", # Light Coral
    "Caudal neurectoderm" : "#FECE44", # Bright Yellow
    "Neural crest" : "#FFC55F", # Yellow-Orange
    "Forebrain/Midbrain/Hindbrain" : "#F89E7B", # Light Coral-Orange
    "Spinal cord" : "#7ECEFD", # Baby Blue
    "Surface ectoderm" : "#C9B1BD", # Pale Mauve
    "Visceral endoderm" : "#E6A0C4", # Light Pink
    "ExE endoderm" : "#E36BAE", # Bright Pink
    "ExE ectoderm" : "#8B5B6E", # Mauve-Brown
    "Parietal endoderm" : "#748CAB", # Blue-Gray
    "Low quality" : "#E5E5E5", # Light Gray
    "Cranial mesoderm" : "#C4C4C4", # Gray
    "Anterior somitic tissues" : "#A4A4A4", # Dark Gray
    "Sclerotome" : "#4D4D4D", # Charcoal
    "Dermomyotome" : "#F8B195", # Dusty Peach
    "Posterior somitic tissues" : "#F67280", # Salmon Pink
    "Presomitic mesoderm" : "#C06C84", # Rose
}

nanostring_cosmx_human_nsclc_cell_type_colors = {
    "tumors" : "#635547",
    "neutrophil" : "#DABE99",
    "T CD8 memory" : "#9e6762",
    "fibroblast" : "#FACB12",
    "B-cell" : "#c19f70",
    "endothelial" : "#0F4A9C",
    "T CD4 memory" : "#F397C0",
    "T CD4 naive" : "#F397C0",
    "NK" : "#EF5A9D",
    "epithelial" : "#EF5A9D",
    "macrophage" : "#C594BF",
    "monocyte" : "#DFCDE4",
    "plasmablast" : "#139992",
    "Treg" : "#3F84AA",
    "T CD8 naive" : "#8DB5CE",
    "mDC" : "#005579",
    "pDC" : "#C9EBFB",
    "mast" : "#C9EBFB"}

nanostring_cosmx_human_nsclc_cell_type_colors = {
    "endothelial": "#FFD8B8", # Peach
    "B-cell": "#F9CB9C", # Light Peach
    "myeloid": "#FFB8A2", # Apricot
    "plasmablasts": "#F08CAE", # Pink
    "fibroblast": "#F49FAD", # Pale Pink
    "neutrophil": "#E7F0A4", # Pale Canary
    "NK/T cell": "#C7F2C2", # Light Green
    "epithelial": "#97E6A1", # Soft Green
    "fibrbolblast": "#FE938C", # Coral
    "mast": "#A5DEE4", # Pale Blue
    "tumors": "#5CA4A9", # Blue-Green
}

vizgen_merfish_mouse_liver_cell_type_colors = {
    "Hepatocyte" : "#635547",
    "Macrophage" : "#DABE99",
    "SEC" : "#9e6762",
    "Erythroid-cell_Erythroid-progenitor_Hepatocyte_MK_Neutrophil" : "#FACB12",
    "AEC_Potential-HSC" : "#c19f70",
    "HSC_Pre-B-cell" : "#0F4A9C",
    "MK" : "#F397C0",
    "AEC_Hepatocyte" : "#F397C0",
    "HSC" : "#EF5A9D",
    "Erythroid-progenitor_Hepatocyte_Neutrophil_Pre-B-cell" : "#EF5A9D",
    "Neutrophil" : "#C594BF",
    "Pre-B-cell" : "#DFCDE4"}

vizgen_merfish_human_ovarian_cancer_cell_type_colors = {
    "Fibroblasts" : "#635547",
    "Epithelial cells" : "#DABE99",
    "Macrophages" : "#9e6762",
    "T cells" : "#FACB12",
    "Endothelial cells" : "#c19f70",
    "B cells" : "#0F4A9C"}

spatial_atac_rna_seq_mouse_embryo_and_brain_rna_colors = {
    "R5" : "#635547",
    "R0" : "#DABE99",
    "R3" : "#9e6762",
    "R1" : "#FACB12",
    "R4" : "#c19f70",
    "R7" : "#0F4A9C",
    "R2" : "#F397C0",
    "R8" : "#F397C0",
    "R12" : "#EF5A9D",
    "R10" : "#EF5A9D",
    "R9" : "#C594BF",
    "R11" : "#DFCDE4",
    "R6" : "#139992",
    "R13" : "#3F84AA"}

spatial_atac_rna_seq_mouse_embryo_and_brain_atac_colors = {
    "A0" : "#635547",
    "A4" : "#DABE99",
    "A5" : "#9e6762",
    "A1" : "#FACB12",
    "A6" : "#c19f70",
    "A7" : "#0F4A9C",
    "A2" : "#F397C0",
    "A3" : "#F397C0",
    "C6" : "#635547",
    "C2" : "#DABE99",
    "C5" : "#9e6762",
    "C0" : "#FACB12",
    "C4" : "#c19f70",
    "C7" : "#0F4A9C",
    "C10" : "#F397C0",
    "C3" : "#A5DEE4",
    "C1" : "#5CA4A9",
    "C9" : "#6B5B95",
    "C8" : "#F7DB6A",
    "C11" : "#D6E4B2",
    "C12" : "#218380",
    "C13" : "#ED553B"}

visium_human_heart_colors = {
    "Adip1" : "#FFD8B8", # Peach
    "Adip2" : "#F9CB9C", # Light Peach
    "Adip3" : "#FFB8A2", # Apricot
    "Adip4" : "#F08CAE", # Pink
    "B" : "#F49FAD", # Pale Pink
    "B_cells" : "#E7F0A4", # Pale Canary
    "B_follicular" : "#C7F2C2", # Light Green
    "B_memory" : "#97E6A1", # Soft Green
    "B_plasma" : "#FE938C", # Coral
    "CD4+T_cytox" : "#A5DEE4", # Pale Blue
    "CD4+T_tem" : "#5CA4A9", # Blue-Green
    "CD4T" : "#6B5B95", # Lavender
    "CD4T_Tfh" : "#F9A03F", # Orange
    "CD4T_Th1" : "#F7DB6A", # Light Yellow
    "CD4T_naive" : "#EEBB4D", # Light Amber
    "CD4T_reg" : "#D6E4B2", # Pale Green
    "CD8+T_cytox" : "#A8DADC", # Pale Cyan
    "CD8+T_tem" : "#3D5A80", # Dark Blue
    "CD8T" : "#3E3F8A", # Navy Blue
    "CD14+Mo" : "#218380", # Teal
    "CD16+Mo" : "#90BE6D", # Soft Green-Yellow
    "DC" : "#FFD369", # Yellow
    "DOCK4+MØ1" : "#ED553B", # Red-Orange
    "DOCK4+MØ2" : "#DA627D", # Mauve
    "EC1_cap" : "#6C5B7B", # Purple
    "EC2_cap" : "#4ECDC4", # Mint
    "EC3_cap" : "#65AADD", # Sky Blue
    "EC4_immune" : "#8FBFE0", # Powder Blue
    "EC5_art" : "#A2D2FF", # Pale Sky Blue
    "EC6_ven" : "#F3C969", # Light Amber-Yellow
    "EC7_atria" : "#EE6C4D", # Light Red-Orange
    "EC8_ln" : "#EC4E20", # Bright Red
    "FB1" : "#D64161", # Dark Pink
    "FB2" : "#FF7A5A", # Bright Coral
    "FB3" : "#E7A977", # Light Coral
    "FB4" : "#FECE44", # Bright Yellow
    "FB5" : "#FFC55F", # Yellow-Orange
    "IL17RA+Mo" : "#F89E7B", # Light Coral-Orange
    "ILC" : "#7ECEFD", # Baby Blue
    "LYVE1+MØ1" : "#C9B1BD", # Pale Mauve
    "LYVE1+MØ2" : "#E6A0C4", # Light Pink
    "LYVE1+MØ3" : "#E36BAE", # Bright Pink
    "MAIT" : "#8B5B6E", # Mauve-Brown
    "Mast" : "#748CAB", # Blue-Gray
    "Meso" : "#E5E5E5", # Light Gray
    "Mo_pi" : "#C4C4C4", # Gray
    "MØ_AgP" : "#A4A4A4", # Dark Gray
    "MØ_mod" : "#4D4D4D", # Charcoal
    "NC1" : "#F8B195", # Dusty Peach
    "NC2" : "#F67280", # Salmon Pink
    "NC3" : "#C06C84", # Rose
    "NC4": "#6C5B7B", # Muted Purple
    "NC5": "#355C7D", # Dark Slate Blue
    "NC6": "#6C7B95", # Gray-Blue
    "NK": "#D6BCC0", # Light Mauve
    "NKT": "#D5B9B2", # Dusty Rose
    "NK_ITGAD": "#A56C7B", # Muted Rose
    "NØ": "#F4A261", # Light Orange
    "PC1_vent": "#F29E4C", # Bright Orange
    "PC2_atria": "#E76F51", # Coral-Red
    "PC3_str": "#DA627D", # Mauve
    "SMC1_basic": "#9B4DCA", # Lavender-Purple
    "SMC2_art": "#8A5B5C", # Brownish Red
    "gdT": "#9A8478", # Beige
    "vCM1": "#FF8C42", # Tangerine
    "vCM2": "#F9A11F", # Bright Orange-Yellow
    "vCM3": "#FACC15", # Golden Yellow
    "vCM4": "#E2E062", # Pale Lime
    "vCM5": "#BADE92", # Soft Lime
}

batch_colors = {
    "batch1": "#FFD8B8", # Peach
    "batch2": "#355C7D", # Dark Slate Blue
    "batch3": "#E76F51", # Coral-Red
    "batch4": "#F08CAE", # Pink
    "batch5": "#F49FAD", # Pale Pink
    "batch6": "#E7F0A4", # Pale Canary
    "batch7": "#C7F2C2", # Light Green
    "batch8": "#97E6A1", # Soft Green
    "batch9": "#FE938C", # Coral
    "batch10": "#A5DEE4", # Pale Blue
    "batch11": "#5CA4A9", # Blue-Green
    "batch12": "#6B5B95", # Lavender
    "batch13": "#F9A03F", # Orange
    "batch14": "#F7DB6A", # Light Yellow
    "batch15": "#EEBB4D", # Light Amber
    "batch16": "#D6E4B2", # Pale Green
    "batch17": "#A8DADC", # Pale Cyan
    "batch18": "#3D5A80", # Dark Blue
    "batch19": "#3E3F8A", # Navy Blue
    "batch20": "#218380", # Teal
    "p22": "#FFD8B8", # Peach
    "embryo1_z2": "#FFD8B8", # Peach
    "embryo1_z5": "#355C7D", # Dark Slate Blue
    "embryo2_z2": "#E76F51", # Coral-Red
    "embryo2_z5": "#F08CAE", # Pink
    "embryo3_z2": "#F49FAD", # Pale Pink
    "embryo3_z5": "#E7F0A4", # Pale Canary
    "patient1": "#FFD8B8", # Peach
    "patient2": "#355C7D", # Dark Slate Blue
    "embryo1": "#FFD8B8", # Peach
    "embryo2": "#355C7D", # Dark Slate Blue
    "embryo3": "#E76F51", # Coral-Red
    "lung5_rep1": "#FFD8B8", # Peach
    "lung5_rep2": "#355C7D", # Dark Slate Blue
    "lung5_rep3": "#E76F51", # Coral-Red
    "lung6": "#F08CAE", # Pink
    "lung9_rep2": "#F49FAD", # Pale Pink
    "lung12": "#E7F0A4", # Pale Canary
    "lung13": "#C7F2C2", # Light Green
}

mapping_entity_colors = {
    "reference" : "#FFD8B8", # Peach
    "query" : "#355C7D", # Dark Slate Blue
}

vizgen_merfish_human_lung_cancer = {
    "Epiblast" : "#FFD8B8", # Peach
    "Primitive Streak" : "#F9CB9C", # Light Peach
    "Caudal epiblast" : "#FFB8A2", # Apricot
    "PGC" : "#F08CAE", # Pink
    "Anterior Primitive Streak" : "#F49FAD", # Pale Pink
    "Notochord" : "#E7F0A4", # Pale Canary
    "Def. endoderm" : "#C7F2C2", # Light Green
    "Definitive endoderm" : "#97E6A1", # Soft Green
    "Gut" : "#FE938C", # Coral
    "Gut tube" : "#A5DEE4", # Pale Blue
    "Nascent mesoderm" : "#5CA4A9", # Blue-Green
    "Mixed mesoderm" : "#6B5B95", # Lavender
    "Intermediate mesoderm" : "#F9A03F", # Orange
    "Caudal Mesoderm" : "#F7DB6A", # Light Yellow
    "Paraxial mesoderm" : "#EEBB4D", # Light Amber
    "Somitic mesoderm" : "#D6E4B2", # Pale Green
    "Pharyngeal mesoderm" : "#A8DADC", # Pale Cyan
    "Splanchnic mesoderm" : "#3D5A80", # Dark Blue
    "Cardiomyocytes" : "#3E3F8A", # Navy Blue
    "Allantois" : "#218380", # Teal
    "ExE mesoderm" : "#90BE6D", # Soft Green-Yellow
    "Lateral plate mesoderm" : "#FFD369", # Yellow
    "Mesenchyme" : "#ED553B", # Red-Orange
    "Mixed mesenchymal mesoderm" : "#DA627D", # Mauve
    "Haematoendothelial progenitors" : "#6C5B7B", # Purple
    "Endothelium" : "#4ECDC4", # Mint
    "Blood progenitors 1" : "#65AADD", # Sky Blue
    "Blood progenitors 2" : "#8FBFE0", # Powder Blue
    "Erythroid1" : "#A2D2FF", # Pale Sky Blue
    "Erythroid2" : "#F3C969", # Light Amber-Yellow
    "Erythroid3" : "#EE6C4D", # Light Red-Orange
    "Erythroid" : "#EC4E20", # Bright Red
    "Blood progenitors" : "#D64161", # Dark Pink
    "NMP" : "#FF7A5A", # Bright Coral
    "Rostral neurectoderm" : "#E7A977", # Light Coral
    "Caudal neurectoderm" : "#FECE44", # Bright Yellow
    "Neural crest" : "#FFC55F", # Yellow-Orange
    "Forebrain/Midbrain/Hindbrain" : "#F89E7B", # Light Coral-Orange
    "Spinal cord" : "#7ECEFD", # Baby Blue
    "Surface ectoderm" : "#C9B1BD", # Pale Mauve
    "Visceral endoderm" : "#E6A0C4", # Light Pink
    "ExE endoderm" : "#E36BAE", # Bright Pink
    "ExE ectoderm" : "#8B5B6E", # Mauve-Brown
    "Parietal endoderm" : "#748CAB", # Blue-Gray
    "Low quality" : "#E5E5E5", # Light Gray
    "Cranial mesoderm" : "#C4C4C4", # Gray
    "Anterior somitic tissues" : "#A4A4A4", # Dark Gray
    "Sclerotome" : "#4D4D4D", # Charcoal
    "Dermomyotome" : "#F8B195", # Dusty Peach
    "Posterior somitic tissues" : "#F67280", # Salmon Pink
    "Presomitic mesoderm" : "#C06C84", # Rose
}

['Multiciliated (non-nasal)',
 'Classical monocytes',
 'AT2',
 'Suprabasal',
 'Adventitial fibroblasts',
 'Pericytes',
 'SMG serous (bronchial)',
 'EC aerocyte capillary',
 'Migratory DCs',
 'CD4 T cells',
 'CD8 T cells',
 'Hillock-like',
 'Lymphatic EC mature',
 'Goblet (nasal)',
 'Basal resting',
 'Mast cells',
 'Monocyte-derived Mph',
 'Peribronchial fibroblasts',
 'Lymphatic EC proliferating',
 'SMG serous (nasal)',
 'EC venous systemic',
 'EC venous pulmonary',
 'Club (nasal)',
 'DC2',
 'Alveolar fibroblasts',
 'Plasma cells',
 'Multiciliated (nasal)',
 'Mesothelium',
 'B cells',
 'NK cells',
 'Goblet (bronchial)',
 'Smooth muscle',
 'Alveolar Mph CCL3+',
 'Alveolar macrophages',
 'EC general capillary',
 'Lymphatic EC differentiating',
 'Club (non-nasal)',
 'Non-classical monocytes',
 'AT2 proliferating',
 'Interstitial Mph perivascular',
 'AT0',
 'Alveolar Mph proliferating',
 'Tuft',
 'Neuroendocrine',
 'Myofibroblasts',
 'SM activated stress response',
 'Subpleural fibroblasts',
 'SMG mucous',
 'Alveolar Mph MT-positive',
 'EC arterial',
 'Plasmacytoid DCs',
 'Smooth muscle FAM83D+',
 'SMG duct',
 'pre-TB secretory',
 'Deuterosomal',
 'AT1',
 'Hematopoietic stem cells',
 'Goblet (subsegmental)',
 'T cells proliferating',
 'Ionocyte',
 'DC1']