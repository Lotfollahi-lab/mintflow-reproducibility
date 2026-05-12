#!/usr/bin/env python
# coding: utf-8

# The main NBscript to run NiCo.

# In[ ]:

PATH_1 = "/nfs/team361/aa36/OnGit/nico_tutorial/NiCo/"
PATH_2 = "/nfs/team361/aa36/OnGit/nico_tutorial/"

import os, sys
sys.path.append(PATH_1)
sys.path.append(PATH_2)


# if you installed the nico package 

import argparse
import types
import yaml
import anndata
import pickle

import NiCo as nico 
from NiCo import Annotations as sann
from NiCo import Interactions as sint
from NiCo import Covariations as scov

import numpy as np
import os
import matplotlib.pyplot as plt 
from matplotlib.collections import PatchCollection
import time
import scanpy as sc
import argparse


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


args = types.SimpleNamespace()
if(flag_isrunning_jupyterNB == False):
    parser = argparse.ArgumentParser(description='dsc')

# ds_name
if(flag_isrunning_jupyterNB):
    args.ds_name = "1024HVGs_2_Psoriasis"
else:
    parser.add_argument('--ds_name', type=str, help = 'ddd.')

# str_prefix_runname
if(flag_isrunning_jupyterNB):
    args.str_prefix_runname = "Attempt_4Debug_Delete_44"
else:
    parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')


# num_factors
if(flag_isrunning_jupyterNB):
    args.num_factors = 5
else:
    parser.add_argument('--num_factors', type=int, help = 'ddd.')


if(flag_isrunning_jupyterNB == False):
    args = parser.parse_args()
print("args = {}".format(args)) #======================================================


# In[ ]:


str_runname = "ds_{}_alltissuescombined_{}".format(
    args.ds_name,
    args.str_prefix_runname
)
str_runname


# In[ ]:


# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# print(device)


# In[ ]:





# # 1. Load the yaml files (ds info, etc)

# In[ ]:


with open("../experiments.yml", 'r') as f:
    content_experiments_yml = yaml.safe_load(f)


# In[ ]:


with open("../config_ds_annotations.yml", 'r') as f:
    content_config_ds_annotations_yml = yaml.safe_load(f)


# In[ ]:





# # 2. Create an anndata object, containning all sections for NiCo

# In[ ]:


content_experiments_yml[args.ds_name]


# In[ ]:


adata = anndata.concat(
    [sc.read_h5ad(fname_adata) for fname_adata in content_experiments_yml[args.ds_name]]
)
adata


# In[ ]:


content_config_ds_annotations_yml[args.ds_name]


# In[ ]:


path_myworkspace = "./NonGit/Runs_TisssuesCombined/{}/".format(str_runname)
os.makedirs(path_myworkspace, exist_ok=True)

path_workspace = path_myworkspace

path_nicoworkspace = path_myworkspace  # os.path.join(path_myworkspace, 'NiCoOutput/')
# os.makedirs(path_nicoworkspace, exist_ok=True)


# In[ ]:


adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
].cat.remove_unused_categories()

# adata.obs_names = adata.obs_names.astype(str)

adata.obs_names = "CELL_" + adata.obs_names.astype(str)  # as if otherwise it's a problem for nico


# In[ ]:


adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
].value_counts()


# In[ ]:


set(
    adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']].tolist()
)


# In[ ]:


# adata = adata[
#     adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']].isin([
#          'output-XETG00155__0056910__CE5-SKI-28-FO-1-S22-A1__20250321__130447'
#     ]).tolist()
# ].copy()  # ---override TODO:revert, TODO:revert
# adata


# In[ ]:


np.allclose(
    adata.X.data,
    np.floor(adata.X.data)
)
adata.layers['counts'] = adata.X.copy()

adata_original = adata.copy()


# In[ ]:


# some preprocessing, so `adata.uns` is not empty.
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=adata.shape[1], subset=True)
sc.tl.pca(adata)
# sc.pp.neighbors(adata)
sc.pp.neighbors(adata, use_rep='X')  # TODO:provide the tissue section ID
adata.raw = adata


# In[ ]:


# remove cell type labels with very low counts
annotation_slot = content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
cell_type_counts = adata.obs[annotation_slot].value_counts()
singletons = cell_type_counts[cell_type_counts <= 10].index
majority_type = adata.obs[annotation_slot].mode()[0]
adata.obs[annotation_slot] = adata.obs[annotation_slot].replace(singletons, majority_type)
adata.obs[annotation_slot] = adata.obs[annotation_slot].astype('category')

# Replaces common illegal Excel characters in the 'cell_type' column
adata.obs[annotation_slot] = adata.obs[annotation_slot].str.replace(r'[:\\/?*\[\]]', '_', regex=True)

adata


# In[ ]:


adata.obs[annotation_slot].value_counts()


# add .obsm['spatial'] ==================
adata.obs['x'] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_X']
]
adata.obs['y'] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_Y']
]

adata.obsm['spatial'] = np.stack(
    [np.array(adata.obs['x'].tolist()), np.array(adata.obs['y'].tolist())],
    0
).T



adata.write_h5ad(
    os.path.join(
        path_myworkspace,
        'adata.h5ad'
    )
)

adata.write_h5ad(
    os.path.join(
        path_nicoworkspace,
        'adata.h5ad'
    )
)

# sct_spatial

adata.write_h5ad(
    os.path.join(
        path_myworkspace,
        'sct_spatial.h5ad'
    )
)

adata.write_h5ad(
    os.path.join(
        path_nicoworkspace,
        'sct_spatial.h5ad'
    )
)

print("Done!")


# In[ ]:





# # 3. Train NiCo

# In[ ]:





# In[ ]:


annotation_slot = content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
showit=True
saveas='png'
transparent_mode = False


# In[ ]:


# parameters of the nico 
inputRadius=0  # ak: recall: Radius=0 means jaxtacrine signalling, 
# annotation_slot='rctd_first_type' #spatial cell type slot, ak: it's the obskey of cell types.
# fname_anndata = ''
do_not_use_following_CT_in_niche=[]


# In[ ]:


niche_pred_output=sint.spatial_neighborhood_analysis(
    Radius=inputRadius,
    output_nico_dir=path_workspace,
    anndata_object_name='adata.h5ad',
    spatial_cluster_tag=annotation_slot,
    removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

t_begin = time.time()

cov_out=scov.gene_covariation_analysis(
    quepath=path_workspace,
    refpath=path_workspace,
    iNMFmode=False,
    Radius=inputRadius,
    no_of_factors=args.num_factors,
    spatial_integration_modality='single',
    anndata_object_name='adata.h5ad',
    output_niche_prediction_dir=path_workspace,
    ref_cluster_tag=annotation_slot,
    LRdbFilename='NiCoLRdb.txt'
)  # according to gAI, this is related to the NMF step to infer latent factors.

print("Took {} seconds.".format(
    time.time() - t_begin
))


# In[ ]:





# # 4. Get Xmic from NiCo

# ## 4.1. Get the latent factors for all cell

# In[ ]:


data1=np.load(cov_out.outputname,allow_pickle=True)
data1=data1['weighted_neighborhood_of_factors_in_niche']
data=np.nan_to_num(data1)

featureVector=range(cov_out.no_of_pc,data.shape[1]) # #just neighborhood
X_latent_factors = data[:,0:cov_out.no_of_pc]  # of shape [num_cells x num_factors]
print(X_latent_factors.shape)


# In[ ]:


dict_ct_to_idxct = {v:k for k, v in niche_pred_output.nameOfCellType.items()}
dict_ct_to_idxct


# ## 4.2. Decode Xmic Per Cell Type

# In[ ]:





# In[ ]:


list_adata_todump = []
for ct in set(adata.obs[annotation_slot].tolist()):
    idx_ct = dict_ct_to_idxct[ct]
    stdnorm_ct, W_ct = cov_out.ADDED_dict_clidsp_to_stdnorm_and_W[str(idx_ct)]  # [num_genes], [num_genes x num_factors]

    # subset the mic part of W_ct
    percent_variance_explained = np.array(cov_out.save_reg_coef[idx_ct][8])  # to be used to split the factors to Xint and Xmic
    W_ct = W_ct[
        :,
        (percent_variance_explained > np.median(percent_variance_explained)).tolist() 
    ]  # [num_genes x num_somefactors]
    if len(W_ct.shape) == 1:
        W_ct = np.expand_dims(W_ct, -1)


    X_latent_factors_ct = X_latent_factors[
        (adata.obs[annotation_slot] == ct).tolist(),
        :
    ]  # [num_cells_CT x num_factors]
    X_latent_factors_ct = X_latent_factors_ct[
        :,
        (percent_variance_explained > np.median(percent_variance_explained)).tolist()
    ]  # [num_cells_CT x num_some_factors]

    if len(X_latent_factors_ct.shape) == 1:
        X_latent_factors_ct = np.expand_dims(X_latent_factors_ct, -1)

    X_mic_ct = X_latent_factors_ct @ (W_ct.T)  # [num_cells_CT x num_genes]

    # un-normalise X_mic_ct
    X_mic_ct = X_mic_ct * np.expand_dims(stdnorm_ct, 0)
    X_mic_ct = np.expm1(X_mic_ct)

    # create adata_ct
    adata_ct = adata_original[
        (adata.obs[annotation_slot] == ct).tolist(),
        :
    ].copy()
    orig_rowsum = np.expand_dims(
        np.array(adata_ct.X.sum(1).tolist()).flatten(),
        -1
    )
    X_mic_ct = X_mic_ct * (orig_rowsum / 10000.0)
    adata_ct.obsm['NiCo_Xmic'] = X_mic_ct
    list_adata_todump.append(adata_ct)


adata_todump = anndata.concat(list_adata_todump)
adata_todump.X = adata_todump.layers['counts'].copy()
adata_todump


# In[ ]:





# In[ ]:


np.array(adata.X.sum(1).tolist()).flatten()


# # 5. Dump the result

# In[ ]:


path_workspace


# In[ ]:


adata_todump.write_h5ad(
    os.path.join(
        path_workspace,
        'adata_NiCo_Result.h5ad'
    )
)
print("Done!")


# In[ ]:


path_workspace


# In[ ]: