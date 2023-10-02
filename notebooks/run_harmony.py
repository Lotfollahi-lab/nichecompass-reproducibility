import pandas as pd
import numpy as np
import harmonypy as hm

meta_data = pd.read_csv('../datasets/srt_data/irene_temp/nanostring_cosmx_human_nsclc/metadata_all.csv')
vars_use = ['sample']
data_mat = np.array(pd.read_csv('../datasets/srt_data/irene_temp/nanostring_cosmx_human_nsclc/pca_scaled_all.csv'))

ho = hm.run_harmony(data_mat, meta_data, vars_use)
res = pd.DataFrame(ho.Z_corr)
res = res.T
res.columns = ['harmony{}'.format(i + 1) for i in range(res.shape[1])]
res.index = meta_data.index
res.to_csv("../datasets/srt_data/irene_temp/nanostring_cosmx_human_nsclc/harmony_scaled_all.csv", index = False)
