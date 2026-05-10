

import os, sys
import argparse
import types
import yaml
# import mintflow
import pickle
import anndata
from tqdm.autonotebook import tqdm
import time
import scanpy as sc
import re


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd
import numpy as np
from scipy import sparse

# =======
parser = argparse.ArgumentParser(description='dsc')
parser.add_argument('--fname_adata', type=str, help = 'ddd.')
parser.add_argument('--fname_adata_int', type=str, help = 'ddd.')
parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')
parser.add_argument('--flag_drop_homogregions', type=str, help = 'ddd.')
args = parser.parse_args()
print("args = {}".format(args)) #======================================================

# check/preprocess args ====
assert isinstance(args.flag_drop_homogregions, str)
assert (args.flag_drop_homogregions in ['True', 'False'])
args.flag_drop_homogregions = (args.flag_drop_homogregions == 'True')


str_runname = "{}_{}_{}".format(
    args.str_prefix_runname,
    args.fname_adata.split('/')[-2],
    str(args.flag_drop_homogregions)
)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


adata = sc.read_h5ad(
    args.fname_adata
)

adata.obsm['spatial'] = np.stack(
    [np.array(adata.obs['x'].tolist()), np.array(adata.obs['y'].tolist())],
    0
).T

adata.obs_names = "CELL_" + adata.obs_names.astype(str)  # as if otherwise it's a problem for nico
assert np.allclose(
    adata.X.data,
    np.floor(adata.X.data)
)

adata.layers['counts'] = adata.X.copy()

adata_original = adata.copy()


adata_unnorm = adata.copy()


adata_int = sc.read_h5ad(
    args.fname_adata_int
)


# discard the homogenous regions if needed
if args.flag_drop_homogregions:
    adata = adata[
        (adata.obs['region_label'] != "Region5") & (adata.obs['region_label'] != "Region6") 
    ]

if args.flag_drop_homogregions:
    adata_int = adata_int[
        (adata_int.obs['region_label'] != "Region5") & (adata_int.obs['region_label'] != "Region6") 
    ]

if args.flag_drop_homogregions:
    adata_original = adata_original[
        (adata_original.obs['region_label'] != "Region5") & (adata_original.obs['region_label'] != "Region6") 
    ]

if args.flag_drop_homogregions:
    adata_unnorm = adata_unnorm[
        (adata_unnorm.obs['region_label'] != "Region5") & (adata_unnorm.obs['region_label'] != "Region6") 
    ]



import numpy as np
from scipy.stats import wasserstein_distance
from scipy import stats

def func_mse(a, b):
    return 'MSE', np.mean((a-b)**2)

def func_mae(a, b):
    return 'MAE', np.mean(np.abs(a-b))

def func_wassdist(a, b):
    return 'EMD', wasserstein_distance(a.flatten(), b.flatten())

def func_pearsoncorrel(a, b):
    try:
        val_corelcoef = stats.pearsonr(a.flatten(), b.flatten()).statistic
    except:
        val_corelcoef = None
    return 'PearsonCorrelation', val_corelcoef



class EvalLargeReadoutsXsplpredExactVersion:
    '''
    Evaluates predXspl on large readouts (i.e. after excluding small readouts) and when the number of readouts is "exactly" equal to mincut_readout.
    '''

    def __init__(self, mincut_readout:int):
        self.mincut_readout = mincut_readout
        self.list_measures = [func_mse, func_mae, func_wassdist, func_pearsoncorrel]

    def eval(self, np_xspl_gt:np.ndarray, np_xspl_pred:np.ndarray, np_xobs:np.ndarray, flag_normalize:bool):
        
        if not isinstance(np_xspl_gt, np.ndarray):
            assert sparse.issparse(np_xspl_gt)
            np_xspl_gt = np_xspl_gt.toarray()

        if not isinstance(np_xspl_pred, np.ndarray):
            assert sparse.issparse(np_xspl_pred)
            np_xspl_pred = np_xspl_pred.toarray()
        
        if not isinstance(np_xobs, np.ndarray):
            assert sparse.issparse(np_xobs)
            np_xobs = np_xobs.toarray()
        

        set_cnts = list(
            set(np_xobs[np_xobs >= self.mincut_readout].flatten().tolist())
        )
        set_cnts.sort()

        dict_toret = {}
        for min_count in set_cnts:
            # mask_min_exp = (np_xobs >= min_count)
            mask_selecteval = (np_xobs == min_count)
            np_pred = np_xspl_pred + 0.0  # np_xspl_pred[mask_nonzero_exp].flatten() + 0.0
            if flag_normalize:
                try:
                    np_pred = np_pred - np.expand_dims(np.min(np_pred, 1), 1)
                    np_pred = np_pred / np.expand_dims(np.max(np_pred, 1), 1)
                    np_pred = np_pred[mask_selecteval].flatten() * np_xobs[mask_selecteval].flatten()
                except:
                    np_pred = np_xspl_pred[mask_selecteval].flatten() + 0.0
            else:
                np_pred = np_xspl_pred[mask_selecteval].flatten() + 0.0

            np_gt = np_xspl_gt[mask_selecteval].flatten() + 0.0

            for measure in self.list_measures:
                measname, measval = measure(np_pred, np_gt)
                dict_toret["{} (among readout == {}, total={})".format(
                    measname, min_count, np.sum(np_xobs >= min_count))
                ] = measval

        return dict_toret
    

mincut_readout = 20
evaluator = EvalLargeReadoutsXsplpredExactVersion(mincut_readout=mincut_readout)



# create Xmic
Xmic_random_binary = adata.X.copy()
Xmic_random_binary.data = np.round(
    Xmic_random_binary.data * ((np.random.rand(Xmic_random_binary.data.shape[0]) > 0.5) + 0.0)
)
Xmic_random_binary.eliminate_zeros()

Xmic_random_binary = Xmic_random_binary.toarray()


np_xspl_gt = adata_unnorm.X.toarray() - adata_int.X
dict_output_e = evaluator.eval(
    np_xspl_gt=np_xspl_gt,
    np_xspl_pred=Xmic_random_binary,
    np_xobs=adata_unnorm.X.toarray(),
    flag_normalize=False  # OLD: note that for simiv this flag is set to True.
)

path_output_files = os.path.join(
    './NonGit/',
    'OuputRuns',
    str_runname
)
os.makedirs(path_output_files, exist_ok=True)

dict_todump = {
    'dict_result':dict_output_e
}
with open(os.path.join(path_output_files, 'vistoken.pkl'), 'wb') as f:
    pickle.dump(dict_todump, f)
print("Dumped the result.")



# convert `dict_output_e` to a dataframe
path_ouptput_files = path_output_files
pattern = r"(\w+)\s+\(among readout == ([\d.]+),\s*total=(\d+)\)"
df_result_finegrained = []
for k in dict_output_e.keys():
    match = re.search(pattern, k)
    assert match
    metric, readout, total = match.groups()
    df_result_finegrained.append([
        'Random_Binary',
        metric, readout, total, dict_output_e[k]
    ])

df_result_finegrained = pd.DataFrame(
    df_result_finegrained,
    columns=[
        'method',
        'metric', 'readout', 'total', 'value'
    ]
)

df_result_finegrained = df_result_finegrained[df_result_finegrained['metric'] != 'PearsonCorrelation']

df_result_finegrained.to_csv(
    os.path.join(
        path_ouptput_files,
        'df_result_finegrained.csv'
    )
)

# save the overal metric
df_result_overall = []
for metric in set(df_result_finegrained['metric'].tolist()):
    np_numbers = df_result_finegrained[
        df_result_finegrained['metric'] == metric
    ]['value'].tolist()
    np_numbers = np.log(np.array(np_numbers) + 1e-10).mean()


    df_result_overall.append([
        args.fname_adata,
        metric,
        np_numbers
    ])


df_result_overall = pd.DataFrame(
    df_result_overall,
    columns=[
        'file',
        'metric',
        'value'
    ]
)
df_result_overall.to_csv(
    os.path.join(
        path_ouptput_files,
        'df_result_overal.csv'
    )
)


print("Script finished successfuly!")
