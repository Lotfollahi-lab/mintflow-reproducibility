


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

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

# load the 4 configs
with open('./config_data_train.yml', 'r') as f:
    config_data_train = yaml.safe_load(f)

with open('./config_data_evaluation.yml', 'r') as f:
    config_data_evaluation = yaml.safe_load(f)

with open('./config_model.yml', 'r') as f:
    config_model = yaml.safe_load(f)

with open('./config_training.yml', 'r') as f:
    config_training = yaml.safe_load(f)


# verify the 4 configs
config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train) 
config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  
config_training = mintflow.verify_and_postprocess_config_training(config_training) 
print("Finished verifying the 4 configuration objects.")


dict_all4_configs = {
    'config_data_train':config_data_train,
    'config_data_evaluation':config_data_evaluation,
    'config_model':config_model,
    'config_training':config_training
}

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
print("Done!")


path_output_files = "./NonGit/OutputPath/"
os.makedirs(path_output_files, exist_ok=True)



for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    '''
    IMPORTANT NOTE: To change the number of epochs, set `config_training['num_training_epochs']` in previous cells of this notebook
    and please refrain from changing the for loop here to, e.g., `for index_epoch in tqdm(range(10), ...)`.
    Because MintFlow's annealing module presumes that the number of epochs equals `config_training['num_training_epochs']`.
    ''' 
    
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

    # evaluate the model and save the evaluation result for the checkpoint
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



