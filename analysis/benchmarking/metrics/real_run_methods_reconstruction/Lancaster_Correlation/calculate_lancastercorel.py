
import os, sys
import numpy as np
import pickle
import pandas as pd
import scanpy as sc
import argparse

import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
import yaml

# part arguments ====
parser = argparse.ArgumentParser(description="DDD")
parser.add_argument('--ds_name', type=str, help = 'ddd.')
args = parser.parse_args()



ds_name = args.ds_name

num_last_epochs_to_avg = 10


dict_dsname_to_obskeycelltype = {
    'Xenium_old_kidney':'level_2_cell_type',
    'Xenium_old_beacon': 'broad_celltypes',
    'Slide-seq_kidney': 'scANVI_pred',
    'VisiumHD_CRC': 'DeconvolutionLabel1'
}

obskey_celltype = dict_dsname_to_obskeycelltype[ds_name]

# load MCC entropies ===
df_MCC_entropies = pd.read_csv(
    os.path.join(
        "./NonGit/PrecomputedMCCEntropies/",
        ds_name,
        obskey_celltype,
        'df.csv'
    )
)
df_MCC_entropies

dict_ct_to_MCCentropy = {
    df_MCC_entropies.iloc[idx_row]['cell_type']: df_MCC_entropies.iloc[idx_row]['MCC_entropy']
    for idx_row in range(df_MCC_entropies.shape[0])
}
dict_ct_to_MCCentropy

df_results = []
for fname_pkl in os.listdir(
    os.path.join(
        "./NonGit/Output_perCT_R2scores/",
        ds_name,
        obskey_celltype
    )
):
    if True:
        print(fname_pkl)
        with open(
            os.path.join(
                "./NonGit/Output_perCT_R2scores/",
                ds_name,
                obskey_celltype,
                fname_pkl
            ),
            'rb'
        ) as f:
            dict_methodname_to_dict_ct_to_listR2score = pickle.load(f)


        for method_name in dict_methodname_to_dict_ct_to_listR2score.keys():
            for ct in dict_methodname_to_dict_ct_to_listR2score[method_name].keys():
                for idx_u, u in enumerate(
                    dict_methodname_to_dict_ct_to_listR2score[method_name][ct]
                ):
                    print(
                        len(dict_methodname_to_dict_ct_to_listR2score[method_name][ct])
                    )
                    df_results.append([
                        ds_name,
                        method_name,
                        ct,
                        idx_u,
                        u,
                        dict_ct_to_MCCentropy[ct]
                    ])



df_results = pd.DataFrame(
    df_results,
    columns=[
        'dataset',
        'method',
        'cell type',
        'epoch',
        'NLL',
        'MCC entropy of cell type'
    ]
)
df_results

df_results = df_results[
    df_results['epoch'] > (df_results['epoch'].max() - num_last_epochs_to_avg) 
]
df_results

from scipy.stats import rankdata, norm, pearsonr

def lancaster_correlation(x, y):
    n = len(x)
    
    # 1. Convert data to ranks
    ranks_x = rankdata(x)
    ranks_y = rankdata(y)
    
    # 2. Transform ranks to normal (probit) scores
    # scores = norm.ppf(ranks / (n + 1))
    tx = norm.ppf(ranks_x / (n + 1))
    ty = norm.ppf(ranks_y / (n + 1))
    
    # 3. Compute r1 = corr(tx, ty)
    r1, _ = pearsonr(tx, ty)
    
    # 4. Compute r2 = corr(tx^2, ty^2)
    r2, _ = pearsonr(tx**2, ty**2)
    
    # 5. Lancaster Correlation is the max of the absolute values
    return max(abs(r1), abs(r2))


dict_methodname_to_lancaster = dict()
for method_name in list(set(df_results['method'].tolist())):
    r = lancaster_correlation(
        x=df_results[df_results['method']==method_name]['NLL'].to_numpy(),
        y=df_results[df_results['method']==method_name]['MCC entropy of cell type'].to_numpy()
    )
    dict_methodname_to_lancaster[method_name] = r
    # assert False

print("Done!")


df_todump = []
for k, v in dict_methodname_to_lancaster.items():
    df_todump.append([
        k,
        v
    ])

df_todump = pd.DataFrame(
    df_todump,
    columns=[
        'method',
        'lancaster_correlation'
    ]
)
df_todump



df_todump.to_csv(
    os.path.join(
        "./NonGit/Output_perCT_R2scores/",
        ds_name,
        obskey_celltype,
        'lancaster_corel.csv'
    )
)

print("Done!")





