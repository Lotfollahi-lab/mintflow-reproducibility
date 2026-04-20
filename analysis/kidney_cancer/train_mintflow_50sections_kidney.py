import os, sys
import sys
print(sys.executable)

import yaml
import mintflow
import pickle
from tqdm.autonotebook import tqdm
import mintflow.interface.base_interface

import time
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

path_anndata ="/nfs/team361/am74/mintflow/revision/full_data_persection/"
# TODO:MODIFY: set to the path where you've put the `.h5ad` file that you downloaded.

import scanpy as sc
import os
from tqdm import tqdm
import scipy.sparse as sp

files = sorted([f for f in os.listdir(path_anndata) if f.endswith(".h5ad")])

for filename in tqdm(files, desc="Processing .h5ad files"):

    adata = sc.read_h5ad(os.path.join(path_anndata, filename))

    # MintFlow's required ID for tissue section = your "batch"
    adata.obs["TissueSectionID_for_MintFlow"] = adata.obs["batch"].astype("category")

    # MintFlow’s required biological batch = your "donor"
    adata.obs["batchID_for_MintFlow"] = adata.obs["donor"].astype("category")


    # save processed file
    # adata.write_h5ad(os.path.join(path_anndata, f"processed_{filename}"))
adata.obs['TissueSectionID_for_MintFlow']
config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=50,
    num_tissue_sections_evaluation=50
)
# list files
files = sorted([f for f in os.listdir(path_anndata) if f.endswith(".h5ad")])

# config_data_train["list_tissue"] = {}

for idx, filename in enumerate(files, start=1):

    key = f"anndata{idx}"
    filepath = os.path.join(path_anndata, filename)

    # config_data_train["list_tissue"][key] = {}

    config_data_train["list_tissue"][key]["file"] = filepath
    # path to this anndata file

    config_data_train["list_tissue"][key]["obskey_cell_type"] = "level_3_cell_type"
    # column containing cell type labels

    config_data_train["list_tissue"][key]["obskey_sliceid_to_checkUnique"] = "donor"
    # unique slice/tissue ID

    config_data_train["list_tissue"][key]["obskey_x"] = "x_centroid"
    config_data_train["list_tissue"][key]["obskey_y"] = "y_centroid"
    # spatial coordinates

    config_data_train["list_tissue"][key]["obskey_biological_batch_key"] = "batch"
    # batch ID

    # dataloader settings
    config_data_train["list_tissue"][key]["config_dataloader_train"]["width_window"]= 1368

    # neighbourhood graph settings
    config_data_train["list_tissue"][key]["config_neighbourhood_graph"] = {
        "n_neighs": 10,
        "set_diag": "False",
        "delaunay": "False"
    }

print("Config entries created for", len(files), "h5ad files.")


# list all .h5ad files in alphabetical order
files = sorted([f for f in os.listdir(path_anndata) if f.endswith(".h5ad")])

# config_data_evaluation["list_tissue"] = {}

for idx, filename in enumerate(files, start=1):

    key = f"anndata{idx}"
    filepath = os.path.join(path_anndata, filename)

    # config_data_evaluation["list_tissue"][key] = {}

    # Path to this anndata file
    config_data_evaluation["list_tissue"][key]["file"] = filepath

    # Column containing cell type labels
    config_data_evaluation["list_tissue"][key]["obskey_cell_type"] = "level_3_cell_type"

    # unique tissue/slice ID
    config_data_evaluation["list_tissue"][key]["obskey_sliceid_to_checkUnique"] = "donor"

    # spatial coordinates
    config_data_evaluation["list_tissue"][key]["obskey_x"] = "x_centroid"
    config_data_evaluation["list_tissue"][key]["obskey_y"] = "y_centroid"

    # batch identifier
    config_data_evaluation["list_tissue"][key]["obskey_biological_batch_key"] = "batch"

    # dataloader settings for evaluation
    config_data_evaluation["list_tissue"][key]["config_dataloader_test"]["width_window"]=1368

    # neighbourhood graph settings
    config_data_evaluation["list_tissue"][key]["config_neighbourhood_graph"] = {
        "n_neighs": 10,
        "set_diag": "False",
        "delaunay": "False"
    }

print("Evaluation config entries created for", len(files), "h5ad files.")

config_model['coef_xbarint2notbatchID_loss'] = 1.0
config_model['coef_xbarspl2notbatchID_loss'] = 1.0
config_model['coef_flowmatchingloss'] = 0.0
config_model['dict_qname_to_scaleandunweighted'] = "impanddisentgl_int#0.1#True&impanddisentgl_spl#0.0#True&varphi_enc_int#0.0#True&varphi_enc_spl#0.0#True&z#0.0#True&sin#0.0#True&sout#0.0#True"

config_training['num_training_epochs'] = 8

# config_training['num_training_epochs'] = 20
# number of training epochs, i.e. the number of times the model sees the dataset during training.

config_training['flag_use_GPU'] = 'True'
# whether GPU is used.

config_training['flag_enable_wandb'] = 'True'
# if set to True, during training different loss terms are logged to wandb.
# It's highly recommended to enable wandb. Please refer to wandb website for more info: `wandb.ai`
config_training['annealing_decoder_XintXspl_coef_max'] = 0.01
config_training['wandb_project_name'] = 'MintFlow'
# wandb project name (ignored if `config_training['flag_enable_wandb']` is set to False)
config_model['coef_loss_CTpredfromZ'] = 100

config_training['wandb_run_name'] = 'Mintflow_8e_h5ad_01_100'
# wandb run name (ignored if `config_training['flag_enable_wandb']` is set to False)
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
import gc
gc.collect()
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
path_output_files = "/nfs/team361/am74/mintflow/revision/train_model_all_data/window_sz_2280_alldatawithbatcheffect/8epochs_save_h5ad/new_hyperparametes/run2/outputs"
# TODO:MODIFY: the path where checkpoints and other files are saved during training.

for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    '''
    IMPORTANT NOTE: To change the number of epochs, set `config_training['num_training_epochs']` in previous cells of this notebook
    and please refrain from changing the for loop here to, e.g., `for index_epoch in tqdm(range(10), ...)`.
    Because MintFlow's annealing module presumes that the number of epochs equals `config_training['num_training_epochs']`.
    ''' 
    
    # train for one epoch
    trainer.train_one_epoch()

    if True:
        mintflow.interface.base_interface.dump_model(
        model,
        os.path.join(
        "/nfs/team361/am74/mintflow/revision/train_model_all_data/window_sz_2280_alldatawithbatcheffect/8epochs_save_h5ad/new_hyperparametes/run2/outputs",
        "model_epoch_{}.pt".format(index_epoch)
                )
        )

        for idx_tissue in range(len(config_data_train)):
            predictions_tissue = mintflow.predict(
                device=device,
                dict_all4_configs=dict_all4_configs,
                data_mintflow=data_mintflow,
                model=model,
                evalulate_on_sections=[idx_tissue],
            )

            adata_current_tissue = (
                data_mintflow['train_list_tissue_section']
                .list_slice[idx_tissue]
                .adata
            )

            for k in predictions_tissue[f'TissueSection {idx_tissue} (zero-based)']:
                adata_current_tissue.obsm[k] = (
                    predictions_tissue[f'TissueSection {idx_tissue} (zero-based)'][k]
                )

            adata_current_tissue.write_h5ad(
                os.path.join(
                    path_output_files,
                    f'adata_epoch_{index_epoch}_tissuesection_{idx_tissue}.h5ad'
                )
            )

            for k in predictions_tissue[f'TissueSection {idx_tissue} (zero-based)']:
                del adata_current_tissue.obsm[k]

            del predictions_tissue
            gc.collect(); gc.collect(); gc.collect()
            print(f"Done for tissue {idx_tissue}")