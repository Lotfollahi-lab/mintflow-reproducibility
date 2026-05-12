#!/usr/bin/env python
# coding: utf-8

# To tune the window size for the tissue sections used in timing/memory requirement.

# In[ ]:


import time

t_begin_script = time.time()




# In[ ]:


import os, sys
import argparse
import yaml
import mintflow
import scanpy as sc
import pickle
from tqdm.autonotebook import tqdm


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd



# args ========
parser = argparse.ArgumentParser(
    prog="DataProcessor",
    description="A simple script to process files.",
    epilog="Use -h for more information."
)

# general
parser.add_argument("--path_workspace", type=str, required=True, help="ddd")
parser.add_argument("--fname_anndata", type=str, required=True, help="ddd")
parser.add_argument("--obskey_cell_type", type=str, required=True, help="ddd")
parser.add_argument("--width_window", type=int, required=True, help="ddd")
parser.add_argument("--dim_sz", type=int, required=True, help="ddd")
parser.add_argument("--wandb_runname", type=str, required=True, help="ddd")



# parser.add_argument("--annotation_slot", type=str, required=True, help="ddd")
args = parser.parse_args()   # args =======

# settings ===
path_anndata = args.fname_anndata  # '../0_Data/NonGit/2_ProcessedData/NumGenes_1000/adata_numcells_1000000.h5ad'
obskey_cell_type = args.obskey_cell_type  # 'level_2_cell_type'


# In[ ]:


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[ ]:





# In[ ]:


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=1,
    num_tissue_sections_evaluation=1
)


# In[ ]:


# tmp
adata = sc.read_h5ad(
    path_anndata
)
adata


# In[ ]:


adata.obs['level_1_cell_type'].value_counts()


# In[ ]:


adata.obs['level_2_cell_type'].value_counts()


# In[ ]:


adata.obs['level_3_cell_type'].value_counts()


# In[ ]:


adata.obs.columns.tolist()


# In[ ]:


adata.obs['donor']


# In[ ]:


print(
    "X in range [{} , {}], width = {}".format(
        adata.obs['x_centroid'].min(),
        adata.obs['x_centroid'].max(),
        adata.obs['x_centroid'].max() - adata.obs['x_centroid'].min()
    )
)

print(
    "Y in range [{} , {}], width = {}".format(
        adata.obs['y_centroid'].min(),
        adata.obs['y_centroid'].max(),
        adata.obs['y_centroid'].max() - adata.obs['y_centroid'].min()
    )
)


# In[ ]:





# In[ ]:


# configure tissue section 1 =========
config_data_train['list_tissue']['anndata1']['file'] = path_anndata
config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = obskey_cell_type
config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'donor'
config_data_train['list_tissue']['anndata1']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata1']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'donor'
config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = args.width_window
config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 5,
    'set_diag': 'False',
    'delaunay': 'False',
}


# In[ ]:


# configure tissue section 1 =======================
config_data_evaluation['list_tissue']['anndata1']['file'] = path_anndata
config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = obskey_cell_type
config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'donor'
config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'donor'
config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = args.width_window
config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 5,
    'set_diag': 'False',
    'delaunay': 'False',
}


# In[ ]:

config_model['dim_sz'] = args.dim_sz




# In[ ]:


config_training['num_training_epochs'] = 1
config_training['flag_use_GPU'] = 'True'
config_training['flag_enable_wandb'] = 'True'
config_training['wandb_project_name'] = 'MintFlow'
config_training['wandb_run_name'] = args.wandb_runname




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


# In[ ]:


data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)


model = mintflow.setup_model(
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow
)

trainer = mintflow.Trainer(
    dict_all4_configs=dict_all4_configs,
    model=model,
    data_mintflow=data_mintflow
)

path_ouptput_files = args.path_workspace
os.makedirs(path_ouptput_files, exist_ok=True)


t_right_before_training_starts = time.time()


for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    
    # train for one epoch
    trainer.train_one_epoch()


t_after_training = time.time()

# get/save the predictions
predictions = mintflow.predict(
    device=device,
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow,
    model=model,
    evalulate_on_sections="all",
)

predictions = predictions['TissueSection 0 (zero-based)']

adata_todump = data_mintflow['train_list_tissue_section'].list_slice[0].adata

for k, v in predictions.items():
    adata_todump.obsm[k] = v


adata_todump.write_h5ad(
    os.path.join(
        path_ouptput_files,
        'adata_output.h5ad'
    )
)

t_script_finished = time.time()



print("TIMING RESULT *** until training starts *** {} *** training took *** {} *** dumping the results took *** {}".format(
    t_right_before_training_starts - t_begin_script,
    t_after_training - t_right_before_training_starts,
    t_script_finished - t_after_training
))


print("Script finished succesfully!")










