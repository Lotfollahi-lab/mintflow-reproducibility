#!/usr/bin/env python
# coding: utf-8

# Runs mintflow on each dataset, with all provided tissue sections of the dataset.
# 
# The input is only (i) `ds_name`, (ii) the prefix of runname, e.g., "Attempt_1" or so, and (iii) number of epochs, and the rest will be read from the `yaml` files.

# In[ ]:


import os, sys
import argparse
import types
import yaml
import mintflow
import pickle
import anndata
from tqdm.autonotebook import tqdm
import time
import scanpy as sc


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd
import numpy as np
from scipy import sparse


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
    args.ds_name = "4_CRC_VisiumHD_Publis"
else:
    parser.add_argument('--ds_name', type=str, help = 'ddd.')

# str_prefix_runname
if(flag_isrunning_jupyterNB):
    args.str_prefix_runname = "Attempt_4Debug_Delete_2"
else:
    parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')

# # num_epochs
# if(flag_isrunning_jupyterNB):
#     args.num_epochs = 20
# else:
#     parser.add_argument('--num_epochs', type=int, help = 'ddd.')


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


# # 1. Load the yaml files

# In[ ]:


with open("../experiments.yml", 'r') as f:
    content_experiments_yml = yaml.safe_load(f)


# In[ ]:


with open("../config_ds_annotations.yml", 'r') as f:
    content_config_ds_annotations_yml = yaml.safe_load(f)


# In[ ]:


with open("./window_sizes.yml", 'r') as f:
    content_window_sizes_yml = yaml.safe_load(f)


# In[ ]:


with open("./additional_lines_toexec.yml", 'r') as f:
    content_additional_lines_toexec_yml = yaml.safe_load(f)


# In[ ]:



# load all anndata objects in a single one
adata = []
for idx_tissue in range(len(content_experiments_yml[args.ds_name])):
    adata.append(
        sc.read_h5ad(
            content_experiments_yml[args.ds_name][idx_tissue]
        )
    )

adata = anndata.concat(adata)
assert np.allclose(
    adata.X.data,
    np.floor(adata.X.data)
)
assert sparse.issparse(adata.X)




path_ouptput_files = "./NonGit/Runs_TisssuesCombined/{}/".format(str_runname)
os.makedirs(path_ouptput_files, exist_ok=True)

# In[ ]:





# # 7. Dump the predictions

# In[ ]:



adata_todump = adata
Xmic_random_uniform = adata.X.copy()
Xmic_random_uniform.data = np.round(
    Xmic_random_uniform.data * ((np.random.rand(Xmic_random_uniform.data.shape[0]) > 0.5) + 0.0)
)
Xmic_random_uniform.eliminate_zeros()
adata_todump.obsm['Xmic_random_uniform'] = Xmic_random_uniform


# In[ ]:


adata_todump.write_h5ad(
    os.path.join(
        path_ouptput_files,
        'adata_result.h5ad'
    )
)
print("Dumped the predictions!")


# In[ ]:





# In[ ]: