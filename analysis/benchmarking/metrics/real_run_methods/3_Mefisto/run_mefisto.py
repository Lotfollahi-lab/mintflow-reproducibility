

#!/usr/bin/env python
# coding: utf-8

# The main script to run mefisto.

# In[ ]:


# import os

# os.environ["CUPY_NVCC_GENERATE_CODE"] = "current"
# os.environ["CCCL_IGNORE_DEPRECATED_CUDA_BELOW_12"] = "1"
# os.environ['NVCCFLAGS'] = "-allow-unsupported-compiler"


# In[ ]:


import numpy as np
import types

if not hasattr(np, 'mat'):
    np.mat = np.asmatrix

import yaml
import os, sys

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

# import torch


# In[ ]:


# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# print(device)


# In[ ]:


# !which python


# In[ ]:


# !python --version


# In[ ]:


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


# In[ ]:


def func_isrunning_jupyterNB() -> bool:
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

flag_isrunning_jupyterNB = func_isrunning_jupyterNB()
print("flag_isrunning_jupyterNB is set to {}".format(flag_isrunning_jupyterNB))


# In[ ]:


# args ========
args = types.SimpleNamespace()

if(flag_isrunning_jupyterNB == False):
    parser = argparse.ArgumentParser(description='dsc')


# ds_name
if(flag_isrunning_jupyterNB):
    args.ds_name = "1024HVGs_1_RCC"
else:
    parser.add_argument('--ds_name', type=str, help = 'ddd.')


# str_prefix_runname
if(flag_isrunning_jupyterNB):
    args.str_prefix_runname = "Attempt_4Debug_1_DELETE"
else:
    parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')

# num_factors_divby2
if(flag_isrunning_jupyterNB):
    args.num_factors_divby2 = 5
else:
    parser.add_argument('--num_factors_divby2', type=int, help = 'ddd.')

# n_inducing
if(flag_isrunning_jupyterNB):
    args.n_inducing = 100
else:
    parser.add_argument('--n_inducing', type=int, help = 'ddd.')


# index_tissue_section
if(flag_isrunning_jupyterNB):
    args.index_tissue_section = 0
else:
    parser.add_argument('--index_tissue_section', type=int, help = 'ddd.')


if(flag_isrunning_jupyterNB == False):
    args = parser.parse_args()
print("args = {}".format(args)) #======================================================


# In[ ]:





# In[ ]:


str_runname = "ds_{}_tissuesection_{}_alltissuescombined_{}".format(
    args.ds_name,
    args.index_tissue_section,
    args.str_prefix_runname
)
str_runname


# In[ ]:





# # 1. Load the yaml files

# In[ ]:


with open("../experiments.yml", 'r') as f:
    content_experiments_yml = yaml.safe_load(f)


# In[ ]:


with open("../config_ds_annotations.yml", 'r') as f:
    content_config_ds_annotations_yml = yaml.safe_load(f)


# In[ ]:





# # 2. Create an anndata object, containning all sections for mefisto

# In[ ]:


content_experiments_yml[args.ds_name]


# In[ ]:





# In[ ]:


adata = sc.read_h5ad(
    content_experiments_yml[args.ds_name][args.index_tissue_section]
)
adata


# In[ ]:


content_config_ds_annotations_yml[args.ds_name]


# In[ ]:


set(
    adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']].tolist()
)


# In[ ]:


# adata = adata[
#     adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']].isin([
#         'CV2-KID-0-FO-1',
#         'CV2-KID-0-FT-1'
#     ]).tolist()
# ]  # ---override TODO:revert, TODO:revert
# adata


# In[ ]:


# at the following the code modifies `adata.X` ==> saving the raw counts here.
assert np.allclose(
    adata.X.data,
    np.floor(
        adata.X.data
    )
)
adata.layers['counts'] = adata.X.copy()


# # 3. Data setup and training

# In[ ]:


adata.obs['x'] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_X']
]
adata.obs['y'] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_Y']
]


# In[ ]:


adata.obsm['spatial'] = np.stack(
    [np.array(adata.obs['x'].tolist()), np.array(adata.obs['y'].tolist())],
    0
).T


# In[ ]:


adata.obsm['spatial'].shape


# In[ ]:


sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)


# In[ ]:





# In[ ]:


# setup t sdsfdfs
ent = entry_point()


# In[ ]:


mefisto_opts = {
    "model_groups": False,
    "warping": False
}

ent.mefisto_opts = mefisto_opts


# In[ ]:


content_config_ds_annotations_yml[args.ds_name]


# In[ ]:


ent.set_data_options(use_float32=True)


# In[ ]:


# Multi-group setup (grabbed from gAI)
# 1. Identify your groups and views
group_names = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
].unique().tolist()
view_name = "gene_expression"

# 2. Initialize lists for data and covariates
# Data structure: [views][groups]
# Covariates structure: [groups] (each is a matrix of n_samples x n_covariates)
data_list = [[]] 
cov_list = []

for g in group_names:
    # Subset AnnData for the current group
    group_adata = adata[adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']] == g].copy()

    # Extract the expression matrix (ensure it's dense for mofapy2)
    matrix = group_adata.X.toarray() if hasattr(group_adata.X, "toarray") else group_adata.X
    data_list[0].append(matrix)

    # Extract spatial coordinates (usually from .obsm['spatial'])
    # Coordinates must be (n_samples, n_dimensions)
    coords = group_adata.obsm['spatial'].astype(np.float64)
    cov_list.append(coords)

# 3. Pass to the entry point
ent.set_data_matrix(data_list, 
                    views_names=[view_name], 
                    groups_names=group_names)

# 4. Pass the covariates
ent.set_covariates(cov_list)


# In[ ]:


# ent.set_data_from_anndata(adata)#, features_subset="highly_variable")
ent.set_model_options(factors=2*args.num_factors_divby2)


# In[ ]:





# In[ ]:


n_inducing = args.n_inducing# 1000
# ent.set_data_options(
#     groups=adata.obs[
#         content_config_ds_annotations_yml[args.ds_name]['obskey_UID']   
#     ].tolist()
# )
# ent.set_covariates([adata.obsm["spatial"]], covariates_names=["x", "y"])

if flag_GPU_available:
    ent.set_train_options(gpu_device=cp_device, gpu_mode=True)

ent.set_smooth_options(sparseGP=True, frac_inducing=n_inducing/adata.n_obs,
                       start_opt=10, opt_freq=10)


# In[ ]:


# as if scipy.empty does not exist anymore, even in very old versions of scipy.
scipy.empty = np.empty


# In[ ]:





# In[ ]:


ent.build()
print("Right before training")


# In[ ]:


t_begin = time.time()
ent.run()
print("Training took {} seconds.".format(time.time() - t_begin))


# In[ ]:


print("Done!")


# In[ ]:





# # 4. Get the predictions

# In[ ]:


ent.model.nodes["Z"].getExpectations()["E"].shape  # interestingly, shape[0] equals adata.shape[0]


# In[ ]:


ent.model.nodes["W"].getExpectations()[0]["E"].shape  
# interestingly, it's not different in the multi-group setting
#  But isn't it that the 0 index will run up to number of sections ???
#  Even with >1 groups, the retval of .getExpectations is of length 1.


# In[ ]:





# In[ ]:


np_factors = ent.model.nodes["Z"].getExpectations()["E"]
np_weights = ent.model.nodes["W"].getExpectations()[0]["E"]
np_scales_touse = ent.model.train_stats['scales']


# In[ ]:


x = np.matmul(
    np_factors,
    np_weights.T
)
x = np.expm1(x)


# In[ ]:


# here Xint and Xmic are calculated using 
# mintflow-private/Analysis/6_InflowVsPrevmethods_on_ImprovedSrtSimData/VisOutputs/1_NMIARI/3_mefisto.ipynb


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


# In[ ]:


x_int.shape


# # 5. Put in anndata and dump

# In[ ]:


path_ouptput_files = "./NonGit/Runs_TisssuesCombined/{}/".format(str_runname)
os.makedirs(path_ouptput_files, exist_ok=True)


# In[ ]:


adata.uns['mefisto_raw_result'] = dict_checkpoint
adata.X = adata.layers['counts']

assert np.allclose(
    adata.X.data,
    np.floor(
        adata.X.data
    )
)

adata.obsm['mefisto_Xint'] = x_int
adata.obsm['mefisto_Xmic'] = x_spl


# In[ ]:


adata.write_h5ad(
    os.path.join(
        path_ouptput_files,
        'adata_result.h5ad'
    )
)
print("Done!")


# In[ ]:





# In[ ]:





# In[ ]:


np.allclose(x, x_spl)


# In[ ]:


listflag_factor_int


# In[ ]:


dict_checkpoint['np_scales_touse']


# In[ ]:


adata.X.sum(1)


# In[ ]:


adata.layers['counts'].data[0:10]


# In[ ]:

