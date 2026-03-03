#!/usr/bin/env python
# coding: utf-8

# # New notebook to run mintflow on cosmx revised by Akbar - lung12
# 

# NB (to be converted to py) to train mintflow on cosmx.
# 
# The inspection of Dec16th.

# In[ ]:


import os, sys
import yaml
import mintflow
import pickle
from tqdm.autonotebook import tqdm


import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd


# In[ ]:


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[ ]:


mintflow


# # 1. load the default config files

# In[ ]:


list_fname_anndata = [
    '/nfs/team361/ms83/data/lung/cosmx_mintflow/raw_data/Lung6+SMI+Flat+data_rawformintflow.h5ad'
]


# In[ ]:


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=len(list_fname_anndata),
    num_tissue_sections_evaluation=len(list_fname_anndata)
)


# # 2. Make config data train

# In[ ]:


config_data_train['list_tissue']['anndata1']['file'] = list_fname_anndata[0]
config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = 'cell_type'
config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'patient'
config_data_train['list_tissue']['anndata1']['obskey_x'] = 'CenterX_global_px'
config_data_train['list_tissue']['anndata1']['obskey_y'] = 'CenterY_global_px'
config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'patient'
config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = 5000
config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# # 3. Make config data evluation

# In[ ]:

config_data_evaluation['list_tissue']['anndata1']['file'] = list_fname_anndata[0]
config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = 'cell_type'
config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'patient'
config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = 'CenterX_global_px'
config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = 'CenterY_global_px'
config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'patient'
config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = 5000
config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# # 4. Customise config model

# In[ ]:

config_model['coef_flowmatchingloss']=0.0

config_model['dict_qname_to_scaleandunweighted'] =  'impanddisentgl_int#0.1#True&impanddisentgl_spl#0.0#True&varphi_enc_int#0.0#True&varphi_enc_spl#0.0#True&z#0.0#True&sin#0.0#True&sout#0.0#True'


# # 5. Customise config training

# In[ ]:

config_training['num_training_epochs'] = 50

config_training['flag_use_GPU'] = 'True'

config_training['flag_enable_wandb'] = 'True'

config_training['wandb_project_name'] = 'mintflow-lung'

config_training['wandb_run_name'] = 'cosmx_lung6_12_Jan'


# # 6. Verify and process config objects

# In[ ]:


config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train) 


# In[ ]:


config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)


# In[ ]:


config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  


# In[ ]:


config_training = mintflow.verify_and_postprocess_config_training(config_training) 


# In[ ]:


print("Finished verifying the 4 configuration objects.")


# # 7. Setup data/model/trainer

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


path_output_files = "/nfs/team361/ms83/data/lung/cosmx_mintflow/Outputs_nb1_lung6/"
os.makedirs(path_output_files, exist_ok=True)
# TODO:MODIFY: the path where checkpoints and other files are saved during training.


# In[ ]:


for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    
    # train for one epoch
    trainer.train_one_epoch()
    
    # get/save the predictions
    predictions = mintflow.predict(
        device=device,
        dict_all4_configs=dict_all4_configs,
        data_mintflow=data_mintflow,
        model=model,
        evalulate_on_sections="all",
    )
    with open(os.path.join(path_output_files, "predictions_epoch_{}.pkl".format(index_epoch)), 'wb') as f:
        pickle.dump(
            predictions,
            f
        )

    # evaluate the model and save the evaluation result for this checkpoint
    df_evaluation_result = mintflow.evaluate_by_known_signalling_genes(
        device=device,
        dict_all4_configs=dict_all4_configs,
        data_mintflow=data_mintflow,
        model=model,
        evalulate_on_sections='all',
        optional_list_colvaltype_toadd=[['training_epoch', index_epoch, 'category']]
    )
    df_evaluation_result.to_pickle(
        os.path.join(
            path_output_files,
            'df_evaluation_result_epoch_{}.pkl'.format(index_epoch)
        )
    )

    # save the checkpoint
    mintflow.dump_checkpoint(
        model=model,
        data_mintflow=data_mintflow,
        dict_all4_configs=dict_all4_configs,
        path_dump=os.path.join(path_output_files, "checkpoint_epoch_{}.pt".format(index_epoch)),
    )  



