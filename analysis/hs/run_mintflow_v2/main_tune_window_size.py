#!/usr/bin/env python
# coding: utf-8

# This notebook adds the hyperparam settings (window size for now) to the DB for a single tissue section.
# 
# If the field already exists, this script doesn't recompute and simply skips.

# In[1]:


import os, sys
import sqlite3
import argparse
import gc
import time
import types
import yaml
import mintflow
import pickle
from tqdm.autonotebook import tqdm
import json

import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd


list_pathstoadd = [
    "/nfs/team361/aa36/OnGit/inflow-reproducibility/",
    "/nfs/team361/aa36/OnGit/inflow-reproducibility/mintflowprivate/"
]

for path in list_pathstoadd:
    if(path not in sys.path):
        sys.path.append(path)


import mintflowprivate.db_hyperparams.utils


# In[ ]:





# In[2]:


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


# In[3]:


# parse arguments ========================================
args = types.SimpleNamespace()
if(flag_isrunning_jupyterNB == False):
    parser = argparse.ArgumentParser(description='dsc')


# fname_DB
if(flag_isrunning_jupyterNB):
    args.fname_DB = '/nfs/team361/aa36/MintFlowDBs/dataset_XXX_windowsizes.db'
else:
    parser.add_argument(
        '--fname_DB',
        type=str,
        help = 'ddd.'
    )



# str_prefix_dataset
if(flag_isrunning_jupyterNB):
    args.str_prefix_dataset = 'OriginalSubmission_Kidney_Dan_Subpopulation_Analysis'
else:
    parser.add_argument(
        '--str_prefix_dataset',
        type=str,
        help = 'ddd.'
    )


# absfname_anndata
if(flag_isrunning_jupyterNB):
    args.absfname_anndata = \
    '/nfs/team361/aa36/OnGit/inflow-reproducibility/Analysis/48_Kidney_Subpopulation_Rebuttal_Robustness_Test/NonGit/Data/inflow_testdata_RCC2_processed_5K.h5ad'
else:
    parser.add_argument(
        '--absfname_anndata',
        type=str,
        help = 'ddd.'
    )

# str_list_windowsize
if(flag_isrunning_jupyterNB):
    args.str_list_windowsize = "[123, 456]"
else:
    parser.add_argument(
        '--str_list_windowsize',
        type=str,
        help = 'ddd.'
    )


# str_list_pygdltrain_numneighbours
if(flag_isrunning_jupyterNB):
    args.str_list_pygdltrain_numneighbours = "[1, 3]"
else:
    parser.add_argument(
        '--str_list_pygdltrain_numneighbours',
        type=str,
        help = 'ddd.'
    )


# config_neighbourhood_graph
if(flag_isrunning_jupyterNB):
    args.config_neighbourhood_graph = "{'n_neighs': 10, 'set_diag': 'False', 'delaunay': 'False'}"

else:
    parser.add_argument(
        '--config_neighbourhood_graph',
        type=str,
        help = 'ddd.'
    )




# obskey_x
if(flag_isrunning_jupyterNB):
    args.obskey_x = 'x_centroid'
else:
    parser.add_argument(
        '--obskey_x',
        type=str,
        help = 'ddd.'
    )

# obskey_y
if(flag_isrunning_jupyterNB):
    args.obskey_y = 'y_centroid'
else:
    parser.add_argument(
        '--obskey_y',
        type=str,
        help = 'ddd.'
    )


# obskey_cell_type
if(flag_isrunning_jupyterNB):
    args.obskey_cell_type = 'level_2_cell_type_formintflow'
else:
    parser.add_argument(
        '--obskey_cell_type',
        type=str,
        help = 'ddd.'
    )


# obskey_sliceid_to_checkUnique
if(flag_isrunning_jupyterNB):
    args.obskey_sliceid_to_checkUnique = 'batch'
else:
    parser.add_argument(
        '--obskey_sliceid_to_checkUnique',
        type=str,
        help = 'ddd.'
    )


# obskey_biological_batch_key
if(flag_isrunning_jupyterNB):
    args.obskey_biological_batch_key = 'batch'
else:
    parser.add_argument(
        '--obskey_biological_batch_key',
        type=str,
        help = 'ddd.'
    )



if(flag_isrunning_jupyterNB == False):
    args = parser.parse_args()
print("args = {}".format(args)) #======================================================


# In[4]:


# other settings ===
fname_DB = args.fname_DB  #'../../NonGit/DBs/DB_windowsize_tuning.db'



device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[7]:


args.config_neighbourhood_graph.replace(" ", '')


# # 1. Show the tissue size (to be used for determining `args.str_list_windowsize` by user)

# In[8]:


adata = sc.read_h5ad(
    args.absfname_anndata
)
print(
    "X in range: [{}  , {}],    width = {}".format(
        adata.obs[args.obskey_x].min(),
        adata.obs[args.obskey_x].max(),
        adata.obs[args.obskey_x].max() - adata.obs[args.obskey_x].min()
    )
)
print(
    "Y in range: [{}  , {}],    width = {}".format(
        adata.obs[args.obskey_y].min(),
        adata.obs[args.obskey_y].max(),
        adata.obs[args.obskey_y].max() - adata.obs[args.obskey_y].min()
    )
)


# In[ ]:





# In[ ]:





# # 2. Compute the stats for each window width and add them to DB (if they do not already exist)

# In[9]:


assert len(set(adata.obs[args.obskey_sliceid_to_checkUnique].tolist())) == 1
UID_tissuesection = adata.obs[args.obskey_sliceid_to_checkUnique].tolist()[0]


# In[10]:


dict_widthwindow_to_dict_result = {}
for width_window in tqdm(eval(args.str_list_windowsize)):
    for pygdltrain_num_neighbours in tqdm(eval(args.str_list_pygdltrain_numneighbours)):

        """
        The part below was commented out, because some cells may have more or less neighrbors than the specified number.
        """
        # if pygdltrain_num_neighbours > eval(args.config_neighbourhood_graph)['n_neighs']:
        #     print("Num neighbours of the neighbour-loader cannot be larger than that of the neighbourhood graph --> skip?")


        # check if the value exists in DB, and if so, skip
        flag_setting_exists = mintflowprivate.db_hyperparams.utils.check_value_exists(
            db_filename=fname_DB,
            args=args,
            int_w=width_window,
            pygdltrain_num_neighbours=pygdltrain_num_neighbours
        )
        if flag_setting_exists:
            print("{}: this setting already exists in the DB ===> skipped.".format(args, width_window))
            continue




        config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
            num_tissue_sections_training=1,
            num_tissue_sections_evaluation=1
        )

        # make config_data_train and config_data_test (the training and the model configs are left untouched)
        config_data_train['list_tissue']['anndata1']['file'] = args.absfname_anndata
        config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = args.obskey_cell_type
        config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = args.obskey_sliceid_to_checkUnique  
        config_data_train['list_tissue']['anndata1']['obskey_x'] = args.obskey_x
        config_data_train['list_tissue']['anndata1']['obskey_y'] = args.obskey_y
        config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = args.obskey_biological_batch_key
        config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = width_window
        config_data_train['list_tissue']['anndata1']['config_dataloader_train']['num_neighbors'] = pygdltrain_num_neighbours
        config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = eval(args.config_neighbourhood_graph)



        config_data_evaluation['list_tissue']['anndata1']['file'] = args.absfname_anndata
        config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = args.obskey_cell_type
        config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = args.obskey_sliceid_to_checkUnique 
        config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = args.obskey_x
        config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = args.obskey_y
        config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = args.obskey_biological_batch_key
        config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = width_window
        config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = eval(args.config_neighbourhood_graph)


        config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train)
        config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
        config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  
        config_training = mintflow.verify_and_postprocess_config_training(config_training) 


        dict_all4_configs = {
            'config_data_train':config_data_train,
            'config_data_evaluation':config_data_evaluation,
            'config_model':config_model,
            'config_training':config_training
        }


        # get the window size, etc.
        data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)


        dict_widthwindow_to_dict_result[width_window] = {
            'maxnumcentralnodes':data_mintflow['dict_slideID_to_maxnumcentralnodes'][UID_tissuesection] + 0,
            'worsecasebatchsize':data_mintflow['dict_slideID_to_worsecasebatchsize'][UID_tissuesection] + 0
        }

        mintflowprivate.db_hyperparams.utils.save_to_DB(
            db_filename=fname_DB,
            args=args,
            int_w=width_window,
            pygdltrain_num_neighbours=pygdltrain_num_neighbours,
            maxnumcentralnodes=str(data_mintflow['dict_slideID_to_maxnumcentralnodes'][UID_tissuesection] + 0),
            worsecasebatchsize=str(data_mintflow['dict_slideID_to_worsecasebatchsize'][UID_tissuesection] + 0),
        )


        # delete and gc.collect
        del config_data_train, config_data_evaluation, config_model, config_training, dict_all4_configs, data_mintflow
        gc.collect(); gc.collect(); gc.collect(); gc.collect(); gc.collect(); gc.collect();

        if not flag_isrunning_jupyterNB:
            plt.clf()
            plt.close('all')
            gc.collect(); gc.collect(); gc.collect(); gc.collect(); gc.collect(); gc.collect();



print("Done!")

