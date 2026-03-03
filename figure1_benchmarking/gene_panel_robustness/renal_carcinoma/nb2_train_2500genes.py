
# Import libraries
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
import yaml

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

# 1. Load the config files

path_old_kidney_configs = '/nfs/users/nfs_m/ms83/mintflow-revision/old_kidney/old_configs/2500_genes/'

with open(os.path.join(path_old_kidney_configs, 'config_data_train_DJ_v2.yml'), 'r') as f:
   config_data_train = yaml.safe_load(f)

with open(os.path.join(path_old_kidney_configs, 'config_data_test_DJ_v2.yml'), 'r') as f:
   config_data_evaluation = yaml.safe_load(f)

with open(os.path.join(path_old_kidney_configs, 'config_training_DJ_v3.yml'), 'r') as f:
   config_training = yaml.safe_load(f)

with open(os.path.join(path_old_kidney_configs, 'config_model_DJ_v2.yml'), 'r') as f:
   config_model = yaml.safe_load(f)

# 2. Make some changes to the config files (e.g. the name of wandb run) 

config_training['wandb_project_name'] = 'mintflow-kidney'
config_training['wandb_run_name'] = 'old_xenium_2500hvg_20_Feb'

# 3. Verify the configs

config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train)
print("Done!")

config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
print("Done!")

config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  
print("Done!")

config_training = mintflow.verify_and_postprocess_config_training(config_training) 
print("Done!")

# 4. Setup the Data/Model/Trainer

dict_all4_configs = {
    'config_data_train':config_data_train,
    'config_data_evaluation':config_data_evaluation,
    'config_model':config_model,
    'config_training':config_training
}

data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
print("Done!")

model = mintflow.setup_model(
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow
)
print("Done!")

trainer = mintflow.Trainer(
    dict_all4_configs=dict_all4_configs,
    model=model,
    data_mintflow=data_mintflow
)
print("Done!")

# 4. Train the model

path_output_files = "/nfs/team361/ms83/data/kidney/xenium_mintflow/Outputs_nb1_2500genes/"
os.makedirs(path_output_files, exist_ok=True)

for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    trainer.train_one_epoch()
    
    mintflow.dump_model(
        model,
        os.path.join(
            path_output_files,
            'model_checkpoint_epoch_{}.pt'.format(index_epoch)
        )
    )

    if index_epoch >= (config_training['num_training_epochs'] - 1 - 2):
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

        # save the checkpoint
        mintflow.dump_checkpoint(
            model=model,
            data_mintflow=data_mintflow,
            dict_all4_configs=dict_all4_configs,
            path_dump=os.path.join(path_output_files, "checkpoint_epoch_{}.pt".format(index_epoch)),
        )