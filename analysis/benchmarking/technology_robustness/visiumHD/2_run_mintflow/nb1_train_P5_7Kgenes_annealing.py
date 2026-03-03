#!/usr/bin/env python
# coding: utf-8

# # New notebook to run mintflow on VisiumHD revised by Akbar - P5
# 

# NB (to be converted to py) to train mintflow on visium HD.
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
    '/nfs/team361/ms83/data/CRC/visiumHD_mintflow/raw_data/adata_sc_P5_rawformintflow_7K_svg.h5ad',
    # '/nfs/team361/ms83/data/CRC/visiumHD_mintflow/raw_data/adata_sc_P2_rawformintflow_7K_svg.h5ad',
    # '/nfs/team361/ms83/data/CRC/visiumHD_mintflow/raw_data/adata_sc_P5_rawformintflow_7K_svg.h5ad'
]


# In[ ]:


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=len(list_fname_anndata),
    num_tissue_sections_evaluation=len(list_fname_anndata)
)


# # 2. Make config data train

# In[ ]:


config_data_train['list_tissue']['anndata1']['file'] = list_fname_anndata[0]
config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = 'DeconvolutionLabel2'
config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'tissue'
config_data_train['list_tissue']['anndata1']['obskey_x'] = 'X'
config_data_train['list_tissue']['anndata1']['obskey_y'] = 'Y'
config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'patient'
config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = 12
config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'True',
}


# In[ ]:


# config_data_train['list_tissue']['anndata2']['file'] = list_fname_anndata[1]
# config_data_train['list_tissue']['anndata2']['obskey_cell_type'] = 'DeconvolutionLabel1'
# config_data_train['list_tissue']['anndata2']['obskey_sliceid_to_checkUnique'] = 'tissue'
# config_data_train['list_tissue']['anndata2']['obskey_x'] = 'X'
# config_data_train['list_tissue']['anndata2']['obskey_y'] = 'Y'
# config_data_train['list_tissue']['anndata2']['obskey_biological_batch_key'] = 'patient'
# config_data_train['list_tissue']['anndata2']['config_dataloader_train']['width_window'] = 5
# config_data_train['list_tissue']['anndata2']['config_neighbourhood_graph'] = {
#     'n_neighs': 10,
#     'set_diag': 'False',
#     'delaunay': 'False',
# }


# In[ ]:


# config_data_train['list_tissue']['anndata3']['file'] = list_fname_anndata[2]
# config_data_train['list_tissue']['anndata3']['obskey_cell_type'] = 'DeconvolutionLabel1'
# config_data_train['list_tissue']['anndata3']['obskey_sliceid_to_checkUnique'] = 'tissue'
# config_data_train['list_tissue']['anndata3']['obskey_x'] = 'X'
# config_data_train['list_tissue']['anndata3']['obskey_y'] = 'Y'
# config_data_train['list_tissue']['anndata3']['obskey_biological_batch_key'] = 'patient'
# config_data_train['list_tissue']['anndata3']['config_dataloader_train']['width_window'] = 5
# config_data_train['list_tissue']['anndata3']['config_neighbourhood_graph'] = {
#     'n_neighs': 10,
#     'set_diag': 'False',
#     'delaunay': 'False',
# }


# In[ ]:


# config_data_train['list_tissue']['anndata4']['file'] = list_fname_anndata[3]
# config_data_train['list_tissue']['anndata4']['obskey_cell_type'] = 'DeconvolutionLabel1'
# config_data_train['list_tissue']['anndata4']['obskey_sliceid_to_checkUnique'] = 'tissue'
# config_data_train['list_tissue']['anndata4']['obskey_x'] = 'X'
# config_data_train['list_tissue']['anndata4']['obskey_y'] = 'Y'
# config_data_train['list_tissue']['anndata4']['obskey_biological_batch_key'] = 'patient'
# config_data_train['list_tissue']['anndata4']['config_dataloader_train']['width_window'] = 5
# config_data_train['list_tissue']['anndata4']['config_neighbourhood_graph'] = {
#     'n_neighs': 10,
#     'set_diag': 'False',
#     'delaunay': 'False',
# }


# # 3. Make config data evluation

# In[ ]:


config_data_evaluation['list_tissue']['anndata1']['file'] = list_fname_anndata[0]
config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = 'DeconvolutionLabel2'
config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'tissue'
config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = 'X'
config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = 'Y'
config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'patient'
config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = 12
config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'True',
}


# In[ ]:


# config_data_evaluation['list_tissue']['anndata2']['file'] = list_fname_anndata[1]
# config_data_evaluation['list_tissue']['anndata2']['obskey_cell_type'] = 'DeconvolutionLabel1'
# config_data_evaluation['list_tissue']['anndata2']['obskey_sliceid_to_checkUnique'] = 'tissue'
# config_data_evaluation['list_tissue']['anndata2']['obskey_x'] = 'X'
# config_data_evaluation['list_tissue']['anndata2']['obskey_y'] = 'Y'
# config_data_evaluation['list_tissue']['anndata2']['obskey_biological_batch_key'] = 'patient'
# config_data_evaluation['list_tissue']['anndata2']['config_dataloader_test']['width_window'] = 5
# config_data_evaluation['list_tissue']['anndata2']['config_neighbourhood_graph'] = {
#     'n_neighs': 10,
#     'set_diag': 'False',
#     'delaunay': 'False',
# }


# In[ ]:


# config_data_evaluation['list_tissue']['anndata3']['file'] = list_fname_anndata[2]
# config_data_evaluation['list_tissue']['anndata3']['obskey_cell_type'] = 'DeconvolutionLabel1'
# config_data_evaluation['list_tissue']['anndata3']['obskey_sliceid_to_checkUnique'] = 'tissue'
# config_data_evaluation['list_tissue']['anndata3']['obskey_x'] = 'X'
# config_data_evaluation['list_tissue']['anndata3']['obskey_y'] = 'Y'
# config_data_evaluation['list_tissue']['anndata3']['obskey_biological_batch_key'] = 'patient'
# config_data_evaluation['list_tissue']['anndata3']['config_dataloader_test']['width_window'] = 5
# config_data_evaluation['list_tissue']['anndata3']['config_neighbourhood_graph'] = {
#     'n_neighs': 10,
#     'set_diag': 'False',
#     'delaunay': 'False',
# }


# In[ ]:


# config_data_evaluation['list_tissue']['anndata4']['file'] = list_fname_anndata[3]
# config_data_evaluation['list_tissue']['anndata4']['obskey_cell_type'] = 'DeconvolutionLabel1'
# config_data_evaluation['list_tissue']['anndata4']['obskey_sliceid_to_checkUnique'] = 'tissue'
# config_data_evaluation['list_tissue']['anndata4']['obskey_x'] = 'X'
# config_data_evaluation['list_tissue']['anndata4']['obskey_y'] = 'Y'
# config_data_evaluation['list_tissue']['anndata4']['obskey_biological_batch_key'] = 'patient'
# config_data_evaluation['list_tissue']['anndata4']['config_dataloader_test']['width_window'] = 5
# config_data_evaluation['list_tissue']['anndata4']['config_neighbourhood_graph'] = {
#     'n_neighs': 10,
#     'set_diag': 'False',
#     'delaunay': 'False',
# }


# In[ ]:





# # 4. Customise config model

# In[ ]:


# config_model['coef_xbarint2notbatchID_loss'] = 1.0
# config_model['coef_xbarspl2notbatchID_loss'] = 1.0


# # 5. Customise config training

# In[ ]:


config_training['num_training_epochs'] = 20

config_training['flag_use_GPU'] = 'True'

config_training['flag_enable_wandb'] = 'True'

config_training['wandb_project_name'] = 'mintflow-CRC'

config_training['wandb_run_name'] = 'visiumHD_P5_7K_12_Jan'

config_training['annealing_decoder_XintXspl_coef_max'] = 0.01


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


path_output_files = "/nfs/team361/ms83/data/CRC/visiumHD_mintflow/Outputs_nb1_P5_annealing/"
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


