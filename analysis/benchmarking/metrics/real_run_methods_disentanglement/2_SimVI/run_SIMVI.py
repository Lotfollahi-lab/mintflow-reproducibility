
#!/usr/bin/env python
# coding: utf-8

# The main script to run SIMVI

# In[ ]:


import os, sys
import numpy as np
from scipy import sparse
import scanpy as sc
# import squidpy as sq
from simvi.model import SimVI
import torch
from torch_geometric.utils.convert import from_scipy_sparse_matrix
import torch_geometric as pyg

import types
import argparse
import matplotlib.pyplot as plt

import yaml
import anndata

import time

# import mintflow
# from mintflow.evaluation.predxspl import EvalLargeReadoutsXsplpredExactVersion


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





# In[ ]:


args = types.SimpleNamespace()
if(flag_isrunning_jupyterNB == False):
    parser = argparse.ArgumentParser(description='dsc')

# ds_name
if(flag_isrunning_jupyterNB):
    args.ds_name = "1_RCC"
else:
    parser.add_argument('--ds_name', type=str, help = 'ddd.')

# str_prefix_runname
if(flag_isrunning_jupyterNB):
    args.str_prefix_runname = "Attempt_4Debug_Delete_3"
else:
    parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')


# batch_size
if(flag_isrunning_jupyterNB):
    args.batch_size = 500
else:
    parser.add_argument('--batch_size', type=int, help = 'ddd.')

# dim_SZ
if(flag_isrunning_jupyterNB):
    args.dim_SZ = 20
else:
    parser.add_argument('--dim_SZ', type=int, help = 'ddd.')


# mae_epochs
if(flag_isrunning_jupyterNB):
    args.mae_epochs = 0
else:
    parser.add_argument('--mae_epochs', type=int, help = 'ddd.')


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


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[ ]:





# # 1. Load the yaml files (ds info, etc)

# In[ ]:


with open("../experiments.yml", 'r') as f:
    content_experiments_yml = yaml.safe_load(f)


# In[ ]:


with open("../config_ds_annotations.yml", 'r') as f:
    content_config_ds_annotations_yml = yaml.safe_load(f)


# In[ ]:





# # 2. Create an anndata object, containning all sections for SIMVI

# In[ ]:


content_experiments_yml[args.ds_name]


# In[ ]:


adata = anndata.concat(
    [sc.read_h5ad(fname_adata) for fname_adata in content_experiments_yml[args.ds_name]]
)
adata


# In[ ]:


# adata = adata[:, 0:100]  # TODO:--override, TODO:revert, TODO:revert
# adata = adata[
#     adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']] == list(
#         set(
#             adata.obs[content_config_ds_annotations_yml[args.ds_name]['obskey_UID']].tolist()
#         )
#     )[0],
#     :
# ]
# adata # TODO:--override, TODO:revert, TODO:revert


# In[ ]:


content_config_ds_annotations_yml[args.ds_name]


# # 3. Data setup and training

# In[ ]:


adata.obs['x'] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_X']
]
adata.obs['y'] = adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_Y']
]


# In[ ]:


SimVI.setup_anndata(
    adata,
    batch_key=content_config_ds_annotations_yml[args.ds_name]['obskey_UID'],
    labels_key=content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
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
    batch_key=content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
)


# In[ ]:


# convert some `adata.obs` columns to category (to avoid one_hot ... error)
for _, v in content_config_ds_annotations_yml[args.ds_name].items():
    adata.obs[v] = adata.obs[v].astype('category')
    adata.obs[v] = adata.obs[v].cat.codes.astype('int64')  #.astype('category')
    # adata.obs[v] = adata.obs[v].astype(np.int64)

# adata.obs['batch'] = adata.obs['batch'].astype('category')


# In[ ]:


adata.obs[
    content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
]


# In[ ]:


model = SimVI(
    adata,kl_weight=1,kl_gatweight=0.01,lam_mi=1000,permutation_rate=0.5,
    n_spatial=args.dim_SZ,
    n_intrinsic=args.dim_SZ
)


# In[ ]:





# In[ ]:


t_begin = time.time()
train_loss, val_loss = model.train(
    edge_index,
    max_epochs=100,
    batch_size=args.batch_size,
    use_gpu=True,
    mae_epochs=args.mae_epochs
)  # as if it doesn't print anything during the `mae_epochs`, and it looks like it's stuck. But it actually isn't stuck.
t_end = time.time()


# In[ ]:


print("Took {} seconds.".format(
    t_end - t_begin
))


# In[ ]:





# # 4. Compute Z, S, and spatial effect

# In[ ]:


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
    cell_type_label=content_config_ds_annotations_yml[args.ds_name]['obskey_celltype'],
    batch_label=content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
)
spat_effect = np.sum(se_list, axis=0)


# In[ ]:





# In[ ]:


# post-processs spatial effect
spat_effect = np.exp(spat_effect) + 0.0 
spat_effect = spat_effect * ((adata.X.toarray() > 0) + 0.0)

spat_effect = np.clip(
    spat_effect,
    a_min=None,
    a_max=adata.X.toarray()
)  # because spat_effect someimtes becomes larger than adata.X !!!


# In[ ]:





# # 5. Dump the result

# In[ ]:


path_ouptput_files = "./NonGit/Runs_TisssuesCombined/{}/".format(str_runname)
os.makedirs(path_ouptput_files, exist_ok=True)


# In[ ]:


adata.obsm['simvi_z'] = z_simvi
adata.obsm['simvi_s'] = s_simvi
adata.obsm['simvi_SE'] = sparse.csr_array(spat_effect)


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


'''
scvi-tools                1.4.1
'''


# In[ ]:





# In[ ]:





# In[ ]: