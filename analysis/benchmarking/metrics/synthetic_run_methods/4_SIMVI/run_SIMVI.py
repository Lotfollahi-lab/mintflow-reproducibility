

import os, sys
import numpy as np
from scipy import sparse
import scanpy as sc
# import squidpy as sq
from simvi.model import SimVI
import re
import pickle
from datetime import datetime
import torch
from torch_geometric.utils.convert import from_scipy_sparse_matrix
import torch_geometric as pyg
import pandas as pd

import types
import argparse
import matplotlib.pyplot as plt

import yaml
import anndata

import time



device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

# input args ===
parser = argparse.ArgumentParser(description='dsc')
parser.add_argument('--fname_simulated_data', type=str, help = 'ddd.')
parser.add_argument('--fname_intrinsic_part_of_simulated_data', type=str, help = 'ddd.')
parser.add_argument('--flag_drop_homogregions', type=str, help = 'ddd.')
parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')
parser.add_argument('--num_training_epochs', type=int, help = 'ddd.')
parser.add_argument('--obskey_celltype', type=str, help = 'ddd.')
parser.add_argument('--neighgraph_num_neighbours', type=int, help = 'ddd.')
parser.add_argument('--num_factors', type=int, help = 'ddd.')
parser.add_argument('--max_epochs', type=int, help = 'ddd.')
parser.add_argument('--batch_size', type=int, help = 'ddd.')
parser.add_argument('--mae_epochs', type=int, help = 'ddd.')
args = parser.parse_args()
# ===================


# check/preprocess args ====
assert isinstance(args.flag_drop_homogregions, str)
assert (args.flag_drop_homogregions in ['True', 'False'])
args.flag_drop_homogregions = (args.flag_drop_homogregions == 'True')
args.fname_adata = args.fname_simulated_data
args.fname_adata_int = args.fname_intrinsic_part_of_simulated_data



str_runname = "SIMVI_{}_{}_{}".format(
    args.str_prefix_runname,
    args.fname_adata.split('/')[-2],
    str(args.flag_drop_homogregions)
)

path_output_files = os.path.join(
    './NonGit/',
    'RunOutputs',
    str_runname
)
os.makedirs(path_output_files, exist_ok=True)


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

adata.obs['batch'] = 'DUMMY'


SimVI.setup_anndata(
    adata,
    batch_key='batch',
    labels_key='cell_type'
)

if 'spatial' not in adata.obsm.keys():
    list_X = adata.obs['x'].tolist()
    list_Y = adata.obs['y'].tolist()
    adata.obsm['spatial'] = np.stack(
        [np.array(list_X), np.array(list_Y)],
        0
    ).T
    assert adata.obsm['spatial'].shape[0] == adata.shape[0]
    assert adata.obsm['spatial'].shape[1] == 2


edge_index = SimVI.extract_edge_index(
    adata,
    n_neighbors=5,
    batch_key='batch'
)


model = SimVI(
    adata,
    n_batch=len(set(adata.obs['batch'])),
    kl_weight=1,
    kl_gatweight=0.01,
    lam_mi=1000,
    permutation_rate=0.5,
    n_spatial=args.num_factors,
    n_intrinsic=args.num_factors
)
print("Created the simvi model.")

print("Started at: {}".format(
    datetime.now()
))
t_tic = time.time()
train_loss, val_loss = model.train(
    edge_index,
    max_epochs=args.max_epochs,
    batch_size=args.batch_size,
    use_gpu=True,
    mae_epochs=args.mae_epochs
)
print("Took {} seconds.".format(time.time() - t_tic))



# get SimVI's Z and S 
z_simvi = model.get_latent_representation(edge_index,representation_kind='intrinsic', give_mean=True)
s_simvi = model.get_latent_representation(edge_index,representation_kind='interaction', give_mean=True)


# In[ ]:


# compute the spaptial effect
adata_se = adata.copy()  # as the name implies, it's the unnormalised adata (i.e. with raw counts in `adata.X`)
adata_se.layers['orig_counts'] = adata.X.copy()
sc.pp.normalize_total(adata_se)
sc.pp.log1p(adata_se)


# In[ ]:


import numpy as np
if not hasattr(np, 'mat'):
    np.mat = np.asmatrix


# In[ ]:


adata_se.obsm['simvi_z'] = z_simvi
adata_se.obsm['simvi_s'] = s_simvi



se_list, r2_zlist, r2_slist, r2_zpvlist, r2_spvlist, S = model.get_se(
    edge_index=edge_index,
    adata=adata_se,
    num_arch =7,
    Kfold=1,
    transformation='none',
    cell_type_label='cell_type',
    batch_label='batch'
)

spat_effect = np.sum(se_list, axis=0)

# post-processs spatial effect
spat_effect = np.exp(spat_effect) + 0.0 
spat_effect = spat_effect * ((adata.X.toarray() > 0) + 0.0)

spat_effect = np.clip(
    spat_effect,
    a_min=None,
    a_max=adata.X.toarray()
)


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


np_xspl_gt = adata_unnorm.X.toarray() - adata_int.X
dict_output_e = evaluator.eval(
    np_xspl_gt=np_xspl_gt,
    np_xspl_pred=spat_effect,
    np_xobs=adata_unnorm.X.toarray(),
    flag_normalize=False  # OLD: note that for simiv this flag is set to True.
)


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
        'SIMVI',
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






