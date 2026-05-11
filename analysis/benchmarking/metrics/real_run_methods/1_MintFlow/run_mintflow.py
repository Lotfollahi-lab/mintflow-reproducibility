
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


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd


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

# num_epochs
if(flag_isrunning_jupyterNB):
    args.num_epochs = 20
else:
    parser.add_argument('--num_epochs', type=int, help = 'ddd.')


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





# In[ ]:





# # 2. Create the 4 configs

# In[ ]:


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=len(content_experiments_yml[args.ds_name]),
    num_tissue_sections_evaluation=len(content_experiments_yml[args.ds_name])
)


# In[ ]:


# configure data train
for idx_tissue in range(len(content_experiments_yml[args.ds_name])):
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['file'] = content_experiments_yml[args.ds_name][idx_tissue]
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['obskey_cell_type'] = \
        content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['obskey_sliceid_to_checkUnique'] = \
        content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['obskey_x'] = content_config_ds_annotations_yml[args.ds_name]['obskey_X']
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['obskey_y'] = content_config_ds_annotations_yml[args.ds_name]['obskey_Y']
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['obskey_biological_batch_key'] = \
        content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['config_dataloader_train']['width_window'] = \
        content_window_sizes_yml[args.ds_name][idx_tissue][content_experiments_yml[args.ds_name][idx_tissue]]
    config_data_train['list_tissue'][f'anndata{idx_tissue+1}']['config_neighbourhood_graph'] = {
        'n_neighs': 5,
        'set_diag': 'False',
        'delaunay': 'False',
    }


# In[ ]:


# configure data evaluation
for idx_tissue in range(len(content_experiments_yml[args.ds_name])):
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['file'] = content_experiments_yml[args.ds_name][idx_tissue]
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['obskey_cell_type'] = \
        content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['obskey_sliceid_to_checkUnique'] = \
        content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['obskey_x'] = content_config_ds_annotations_yml[args.ds_name]['obskey_X']
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['obskey_y'] = content_config_ds_annotations_yml[args.ds_name]['obskey_Y']
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['obskey_biological_batch_key'] = \
        content_config_ds_annotations_yml[args.ds_name]['obskey_UID']
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['config_dataloader_test']['width_window'] = \
        content_window_sizes_yml[args.ds_name][idx_tissue][content_experiments_yml[args.ds_name][idx_tissue]]
    config_data_evaluation['list_tissue'][f'anndata{idx_tissue+1}']['config_neighbourhood_graph'] = {
        'n_neighs': 5,
        'set_diag': 'False',
        'delaunay': 'False',
    }


# In[ ]:


config_training['num_training_epochs'] = args.num_epochs
config_training['flag_use_GPU'] = 'True'
config_training['flag_enable_wandb'] = 'True'
config_training['wandb_project_name'] = 'MintFlow'
config_training['wandb_run_name'] = 'MintFlow_{}'.format(str_runname)


# In[ ]:





# # 3. Additional modifications to configs (if any)

# In[ ]:


if args.ds_name in content_additional_lines_toexec_yml.keys():
    exec(content_additional_lines_toexec_yml[args.ds_name])
else:
    print("No additional dataset-specific modifiactions to exec.")


# In[ ]:





# # 4. Very the 4 configs

# In[ ]:


config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train) 


# In[ ]:


config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)


# In[ ]:


config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  


# In[ ]:


config_training = mintflow.verify_and_postprocess_config_training(config_training) 


# In[ ]:





# # 5. Create data/model/trainer

# In[ ]:


dict_all4_configs = {
    'config_data_train':config_data_train,
    'config_data_evaluation':config_data_evaluation,
    'config_model':config_model,
    'config_training':config_training
}


# In[ ]:


data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)


# In[ ]:


model = mintflow.setup_model(
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow
)


# In[ ]:


trainer = mintflow.Trainer(
    dict_all4_configs=dict_all4_configs,
    model=model,
    data_mintflow=data_mintflow
)


# In[ ]:


str_runname


# # 6. The actual training

# In[ ]:


path_ouptput_files = "./NonGit/Runs_TisssuesCombined/{}/".format(str_runname)
os.makedirs(path_ouptput_files, exist_ok=True)


# In[ ]:


for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    t_begin = time.time()

    trainer.train_one_epoch()

    t_end = time.time()

    print(">>>***>>>***>>> Training in epoch {} took {} seconds.".format(
        index_epoch,
        t_end - t_begin
    ))


# In[ ]:





# # 7. Dump the predictions

# In[ ]:


predictions = mintflow.predict(
    device=device,
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow,
    model=model,
    evalulate_on_sections="all",
)


# In[ ]:





# In[ ]:


list_adata = []
for idx_sl, sl in enumerate(data_mintflow['train_list_tissue_section'].list_slice):
    adata_curr = sl.adata_before_scppnormalize_total.copy()
    for k, v in predictions[f'TissueSection {idx_sl} (zero-based)'].items():
        adata_curr.obsm[k] = v

    list_adata.append(adata_curr)

adata_todump = anndata.concat(list_adata)
adata_todump


# In[ ]:


adata_todump.write_h5ad(
    os.path.join(
        path_ouptput_files,
        'adata_resul.h5ad'
    )
)
print("Dumped the predictions!")


# In[ ]:





# In[ ]: