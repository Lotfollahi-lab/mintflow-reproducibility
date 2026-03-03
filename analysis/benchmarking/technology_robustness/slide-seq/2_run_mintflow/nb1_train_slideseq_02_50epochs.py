#!/usr/bin/env python
# coding: utf-8
# %%

# # New notebook to run mintflow on slide-seq revised by Akbar - Kidney slide-seq
# 

# NB (to be converted to py) to train mintflow on slide-seq.
# 
# The inspection of Dec16th.

# %%


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


# %%


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# %%


mintflow


# # 1. load the default config files

# %%


list_fname_anndata = [
    '/nfs/team361/ms83/data/kidney/slideseq_mintflow/raw_data/adata_puck_02_rawformintflow.h5ad',
    # '/nfs/team361/ms83/data/CRC/visiumHD_mintflow/raw_data/adata_sc_P2_rawformintflow_7K_svg.h5ad',
    # '/nfs/team361/ms83/data/CRC/visiumHD_mintflow/raw_data/adata_sc_P5_rawformintflow_7K_svg.h5ad'
]


# %%


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=len(list_fname_anndata),
    num_tissue_sections_evaluation=len(list_fname_anndata)
)


# # 2. Make config data train

# %%


config_data_train['list_tissue']['anndata1']['file'] = list_fname_anndata[0]
config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = 'scANVI_pred'
config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'batch'
config_data_train['list_tissue']['anndata1']['obskey_x'] = 'x'
config_data_train['list_tissue']['anndata1']['obskey_y'] = 'y'
config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'batch'
config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = 800
config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# # 3. Make config data evluation

# %%


config_data_evaluation['list_tissue']['anndata1']['file'] = list_fname_anndata[0]
config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = 'scANVI_pred'
config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'batch'
config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = 'x'
config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = 'y'
config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'batch'
config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = 800
config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# # 4. Customise config model
  


# # 5. Customise config training

# %%


config_training['num_training_epochs'] = 50

config_training['flag_use_GPU'] = 'True'

config_training['flag_enable_wandb'] = 'True'

config_training['wandb_project_name'] = 'mintflow-slideseq'

config_training['wandb_run_name'] = 'puck_02_13K_12_Jan'

config_training['annealing_decoder_XintXspl_coef_min'] = 0.00001
config_training['annealing_decoder_XintXspl_coef_max'] = 0.01


# # 6. Verify and process config objects

# %%


config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train) 


# %%


config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)


# %%


config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  


# %%


config_training = mintflow.verify_and_postprocess_config_training(config_training) 


# %%


print("Finished verifying the 4 configuration objects.")


# # 7. Setup data/model/trainer

# %%


dict_all4_configs = {
    'config_data_train':config_data_train,
    'config_data_evaluation':config_data_evaluation,
    'config_model':config_model,
    'config_training':config_training
}


# %%


data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)


# %%


model = mintflow.setup_model(
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow
)


# %%


trainer = mintflow.Trainer(
    dict_all4_configs=dict_all4_configs,
    model=model,
    data_mintflow=data_mintflow
)


# %%


path_output_files = "/nfs/team361/ms83/data/kidney/slideseq_mintflow/Outputs_nb1_50epochs/puck_02/"
os.makedirs(path_output_files, exist_ok=True)
# TODO:MODIFY: the path where checkpoints and other files are saved during training.


# %%


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



