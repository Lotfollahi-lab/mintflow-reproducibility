

import numpy as np
import types

if not hasattr(np, 'mat'):
    np.mat = np.asmatrix

import yaml
import os, sys
import re

# import matplotlib.pyplot as plt 
# from matplotlib.collections import PatchCollection

import time
import scanpy as sc
import argparse
from datetime import datetime
from mofapy2.run.entry_point import entry_point
import pickle
import scipy

import anndata


import cupy
flag_GPU_available = None
try:
    # Check if a GPU device is detected
    device_count = cupy.cuda.runtime.getDeviceCount()
    print(f"GPUs detected: {device_count}")

    # Check compute capability of the first device
    if device_count > 0:
        cp_device = cupy.cuda.Device(0)
        print(f"Device 0: {cp_device}")
        print(f"Compute Capability: {cp_device.compute_capability}")
        flag_GPU_available = True
except Exception as e:
    print(f"GPU not detected or CuPy error: {e}")
    flag_GPU_available = False

assert flag_GPU_available is not None


# input args ===
parser = argparse.ArgumentParser(description='dsc')
parser.add_argument('--fname_simulated_data', type=str, help = 'ddd.')
parser.add_argument('--fname_intrinsic_part_of_simulated_data', type=str, help = 'ddd.')
parser.add_argument('--flag_drop_homogregions', type=str, help = 'ddd.')
parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')
parser.add_argument('--num_factors_divby2', type=int, help = 'ddd.')
parser.add_argument('--n_inducing', type=int, help = 'ddd.')

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


assert np.allclose(
    adata.X.data,
    np.floor(
        adata.X.data
    )
)
adata.layers['counts'] = adata.X.copy()

adata.obsm['spatial'] = np.stack(
    [np.array(adata.obs['x'].tolist()), np.array(adata.obs['y'].tolist())],
    0
).T

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

ent = entry_point()

mefisto_opts = {
    "model_groups": False,
    "warping": False
}
ent.mefisto_opts = mefisto_opts

ent.set_data_options(use_float32=True)

if False:
    group_names = 'DUMMY'
    view_name = "gene_expression"

    data_list = [[]] 
    cov_list = []

    for g in group_names:
        # Subset AnnData for the current group
        group_adata = adata[adata.obs['batch'] == g].copy()

        # Extract the expression matrix (ensure it's dense for mofapy2)
        matrix = group_adata.X.toarray() if hasattr(group_adata.X, "toarray") else group_adata.X
        data_list[0].append(matrix)

        # Extract spatial coordinates (usually from .obsm['spatial'])
        # Coordinates must be (n_samples, n_dimensions)
        coords = group_adata.obsm['spatial'].astype(np.float64)
        cov_list.append(coords)

    breakpoint()
    # 3. Pass to the entry point
    ent.set_data_matrix(data_list, 
                        views_names=[view_name], 
                        groups_names=[group_names])

    # 4. Pass the covariates
    ent.set_covariates(cov_list)


    ent.set_model_options(factors=2*args.num_factors_divby2)

    n_inducing = args.n_inducing# 1000
    # ent.set_data_options(
    #     groups=adata.obs[
    #         content_config_ds_annotations_yml[args.ds_name]['obskey_UID']   
    #     ].tolist()
    # )
    # ent.set_covariates([adata.obsm["spatial"]], covariates_names=["x", "y"])

ent.set_data_from_anndata(adata)
ent.set_covariates([adata.obsm["spatial"]], covariates_names=["x", "y"])

ent.set_model_options(factors=2*args.num_factors_divby2)

if flag_GPU_available:
    ent.set_train_options(gpu_device=cp_device, gpu_mode=True)

ent.set_smooth_options(sparseGP=True, frac_inducing=args.n_inducing/adata.n_obs,
                       start_opt=10, opt_freq=10)

scipy.empty = np.empty

ent.build()
print("Right before training")


breakpoint()

t_begin = time.time()
ent.run()
print("Training took {} seconds.".format(time.time() - t_begin))

np_factors = ent.model.nodes["Z"].getExpectations()["E"]
np_weights = ent.model.nodes["W"].getExpectations()[0]["E"]
np_scales_touse = ent.model.train_stats['scales']


# In[ ]:

dict_checkpoint = {
    'np_factors':np_factors,
    'np_weights':np_weights,
    'np_scales_touse':np_scales_touse
}  # only to match the notebook mentioned above.


# In[ ]:


# create x_spl and x_int ====
listflag_factor_spatial = (dict_checkpoint['np_scales_touse'] >= np.median(dict_checkpoint['np_scales_touse'])).tolist()
listflag_factor_int = (dict_checkpoint['np_scales_touse'] < np.median(dict_checkpoint['np_scales_touse'])).tolist()
x_spl = np.expm1(
    np.matmul(
        dict_checkpoint['np_factors'][:, listflag_factor_spatial],
        dict_checkpoint['np_weights'][:, listflag_factor_spatial].T
    )
)
x_int = np.expm1(
    np.matmul(
        dict_checkpoint['np_factors'][:, listflag_factor_int],
        dict_checkpoint['np_weights'][:, listflag_factor_int].T
    )
)

# TODO:HERE compute Xmic

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
    np_xspl_pred=x_spl,
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






