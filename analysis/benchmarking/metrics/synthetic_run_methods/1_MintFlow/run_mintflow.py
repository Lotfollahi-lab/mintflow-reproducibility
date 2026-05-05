

"""
The old synth run was different from the default configs on May2nd, 2026. 
So this script is similar to old synth run.
"""

import os, sys
import yaml
import mintflow
import pickle
import scanpy as sc
from tqdm.autonotebook import tqdm
import gc
import re


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import pandas as pd

import argparse
import types
from pathlib import Path


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

# input args ===
parser = argparse.ArgumentParser(description='dsc')
parser.add_argument('--fname_simulated_data', type=str, help = 'ddd.')
parser.add_argument('--fname_intrinsic_part_of_simulated_data', type=str, help = 'ddd.')
parser.add_argument('--str_prefix_runname', type=str, help = 'ddd.')
parser.add_argument('--num_training_epochs', type=int, help = 'ddd.')
args = parser.parse_args()
# ===================


# modify the anndata (to make it MintFlow-friendly) ===
path_workspace = os.path.join(
    './NonGit/',
    'Workspaces/',
    'Workspace_{}'.format(
        args.str_prefix_runname
    )
)
os.makedirs(
    path_workspace,
    exist_ok=True
)
adata = sc.read_h5ad(args.fname_simulated_data)
adata.obs['batch'] = 'Slide1'

adata.write_h5ad(
    os.path.join(
        path_workspace,
        'adata.h5ad'
    )
)

del adata
gc.collect(); gc.collect(); gc.collect()


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=1,
    num_tissue_sections_evaluation=1
)


# configure tissue section 1 =========
config_data_train['list_tissue']['anndata1']['file'] = os.path.join(
    path_workspace,
    'adata.h5ad'
)
config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = 'cell_type'
config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'batch'
config_data_train['list_tissue']['anndata1']['obskey_x'] = 'x'
config_data_train['list_tissue']['anndata1']['obskey_y'] = 'y'
config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'batch'
config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = 500
config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 5,
    'set_diag': 'False',
    'delaunay': 'False',
}


config_data_evaluation['list_tissue']['anndata1']['file'] = os.path.join(
    path_workspace,
    'adata.h5ad'
)
config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = 'cell_type'
config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'batch'
config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = 'x'
config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = 'y'
config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'batch'
config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = 500
config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 5,
    'set_diag': 'False',
    'delaunay': 'False',
}



config_training['num_training_epochs'] = args.num_training_epochs
config_training['flag_use_GPU'] = 'True'
config_training['flag_enable_wandb'] = 'True'
config_training['wandb_project_name'] = 'MintFlow'
config_training['wandb_run_name'] = "{}".format(
    args.str_prefix_runname
)


# breakpoint()

# hard-coded modifications to configs
config_model['dim_sz'] = 20
config_training['annealing_decoder_XintXspl_coef_min'] = 1.0
config_training['annealing_decoder_XintXspl_coef_max'] = 1.0


# verify the 4 config objects 
config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train) 
config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  
config_training = mintflow.verify_and_postprocess_config_training(config_training) 

print("The configs were verified")
# breakpoint()

# setup data/model/trainer
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


path_ouptput_files = "./NonGit/OutputRuns/Output_{}".format(
    args.str_prefix_runname,
)
# Create the directory if it doesn't exist
os.makedirs(path_ouptput_files, exist_ok=True)



# predictions = mintflow.predict(
#     device=device,
#     dict_all4_configs=dict_all4_configs,
#     data_mintflow=data_mintflow,
#     model=model,
#     evalulate_on_sections="all",
# )
# predictions['TissueSection 0 (zero-based)']['MintFlow_Xmic (before_sc_pp_normalize_total)']
#  is a sparse matrix of shape [num_cells x num_genes]


import numpy as np
from scipy.stats import wasserstein_distance
from scipy import stats

def func_mse(a, b):
    return 'MSE', np.mean((a-b)**2)

def func_mae(a, b):
    return 'MAE', np.mean(np.abs(a-b))

def func_wassdist(a, b):
    return 'EMD', wasserstein_distance(a.flatten(), b.flatten())

def func_pearsoncorrel(a, b):
    try:
        val_corelcoef = stats.pearsonr(a.flatten(), b.flatten()).statistic
    except:
        val_corelcoef = None
    return 'PearsonCorrelation', val_corelcoef



class EvalLargeReadoutsXsplpredExactVersion:
    '''
    Evaluates predXspl on large readouts (i.e. after excluding small readouts) and when the number of readouts is "exactly" equal to mincut_readout.
    '''

    def __init__(self, mincut_readout:int):
        self.mincut_readout = mincut_readout
        self.list_measures = [func_mse, func_mae, func_wassdist, func_pearsoncorrel]

    def eval(self, np_xspl_gt:np.ndarray, np_xspl_pred:np.ndarray, np_xobs:np.ndarray, flag_normalize:bool):
        
        if not isinstance(np_xspl_gt, np.ndarray):
            assert sparse.issparse(np_xspl_gt)
            np_xspl_gt = np_xspl_gt.toarray()

        if not isinstance(np_xspl_pred, np.ndarray):
            assert sparse.issparse(np_xspl_pred)
            np_xspl_pred = np_xspl_pred.toarray()
        
        if not isinstance(np_xobs, np.ndarray):
            assert sparse.issparse(np_xobs)
            np_xobs = np_xobs.toarray()
        

        set_cnts = list(
            set(np_xobs[np_xobs >= self.mincut_readout].flatten().tolist())
        )
        set_cnts.sort()

        dict_toret = {}
        for min_count in set_cnts:
            # mask_min_exp = (np_xobs >= min_count)
            mask_selecteval = (np_xobs == min_count)
            np_pred = np_xspl_pred + 0.0  # np_xspl_pred[mask_nonzero_exp].flatten() + 0.0
            if flag_normalize:
                try:
                    np_pred = np_pred - np.expand_dims(np.min(np_pred, 1), 1)
                    np_pred = np_pred / np.expand_dims(np.max(np_pred, 1), 1)
                    np_pred = np_pred[mask_selecteval].flatten() * np_xobs[mask_selecteval].flatten()
                except:
                    np_pred = np_xspl_pred[mask_selecteval].flatten() + 0.0
            else:
                np_pred = np_xspl_pred[mask_selecteval].flatten() + 0.0

            np_gt = np_xspl_gt[mask_selecteval].flatten() + 0.0

            for measure in self.list_measures:
                measname, measval = measure(np_pred, np_gt)
                dict_toret["{} (among readout == {}, total={})".format(
                    measname, min_count, np.sum(np_xobs >= min_count))
                ] = measval

        return dict_toret
    

mincut_readout = 20
evaluator = EvalLargeReadoutsXsplpredExactVersion(mincut_readout=mincut_readout)

adata_unnorm = sc.read_h5ad(
    args.fname_simulated_data
)
adata_int = sc.read_h5ad(
    args.fname_intrinsic_part_of_simulated_data
)

np_xspl_gt = adata_unnorm.X.toarray() - adata_int.X



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

Xmic_pred = predictions['TissueSection 0 (zero-based)']['MintFlow_Xmic (before_sc_pp_normalize_total)']

dict_output_e = evaluator.eval(
    np_xspl_gt=np_xspl_gt,
    np_xspl_pred=Xmic_pred.toarray(),
    np_xobs=adata_unnorm.X.toarray(),
    flag_normalize=False  # OLD: note that for simiv this flag is set to True.
)


# breakpoint()

dict_todump = {
    'dict_result':dict_output_e
}
with open(os.path.join(path_ouptput_files, 'vistoken.pkl'), 'wb') as f:
    pickle.dump(dict_todump, f)
print("Dumped the result.")

# convert `dict_output_e` to a dataframe
pattern = r"(\w+)\s+\(among readout == ([\d.]+),\s*total=(\d+)\)"
df_result_finegrained = []
for k in dict_output_e.keys():
    match = re.search(pattern, k)
    assert match
    metric, readout, total = match.groups()
    df_result_finegrained.append([
        'MintFlow',
        metric, readout, total, dict_output_e[k]
    ])

df_result_finegrained = pd.DataFrame(
    df_result_finegrained,
    columns=[
        'method',
        'metric', 'readout', 'total', 'value'
    ]
)


# breakpoint()

df_result_finegrained = df_result_finegrained[df_result_finegrained['metric'] != 'PearsonCorrelation']

df_result_finegrained.to_csv(
    os.path.join(
        path_ouptput_files,
        'df_result_finegrained.csv'
    )
)

# breakpoint()

# save the overal metric
df_result_overall = []
for metric in set(df_result_finegrained['metric'].tolist()):
    np_numbers = df_result_finegrained[
        df_result_finegrained['metric'] == metric
    ]['value'].tolist()
    np_numbers = np.log(np.array(np_numbers) + 1e-10).mean()


    df_result_overall.append([
        args.fname_simulated_data,
        metric,
        np_numbers
    ])

df_result_overall = pd.DataFrame(
    df_result_overall,
    columns=[
        'file',
        'metric',
        'value'
    ]
)
df_result_overall.to_csv(
    os.path.join(
        path_ouptput_files,
        'df_result_overal.csv'
    )
)


print("Script finished successfuly!")
