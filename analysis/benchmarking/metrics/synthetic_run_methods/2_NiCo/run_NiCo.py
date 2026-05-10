



import os, sys
PATH_1 = "/nfs/team361/aa36/OnGit/nico_tutorial/NiCo/"
PATH_2 = "/nfs/team361/aa36/OnGit/nico_tutorial/"
sys.path.append(PATH_1)
sys.path.append(PATH_2)


# if you installed the nico package 

import argparse
import types
import yaml
import anndata
import pickle
from scipy import sparse
import pandas as pd
import re

import NiCo as nico 
from NiCo import Annotations as sann
from NiCo import Interactions as sint
from NiCo import Covariations as scov

import numpy as np
import os
import matplotlib.pyplot as plt 
from matplotlib.collections import PatchCollection
import time
import scanpy as sc
import argparse

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


args = types.SimpleNamespace()
if(flag_isrunning_jupyterNB == False):
    parser = argparse.ArgumentParser(description='dsc')

# fname_adata
if(flag_isrunning_jupyterNB):
    args.fname_adata = 0
else:
    parser.add_argument('--fname_adata', type=str, help = 'ddd.')


# fname_adata_int
if(flag_isrunning_jupyterNB):
    args.fname_adata_int = 0
else:
    parser.add_argument('--fname_adata_int', type=str, help = 'ddd.')


# flag_drop_homogregions
if(flag_isrunning_jupyterNB):
    args.flag_drop_homogregions = 'True'
else:
    parser.add_argument('--flag_drop_homogregions', type=str,
                        help = 'ddd.')

# num_factors
if(flag_isrunning_jupyterNB):
    args.num_factors = 'True'
else:
    parser.add_argument('--num_factors', type=int,
                        help = 'ddd.')

if(flag_isrunning_jupyterNB == False):
    args = parser.parse_args()
print("args = {}".format(args)) #======================================================


# check/preprocess args ====
assert isinstance(args.flag_drop_homogregions, str)
assert (args.flag_drop_homogregions in ['True', 'False'])
args.flag_drop_homogregions = (args.flag_drop_homogregions == 'True')

str_runname = "NiCo_{}".format(
    args.fname_adata.split('/')[-2]
)
str_runname



adata = sc.read_h5ad(
    args.fname_adata
)

adata.obsm['spatial'] = np.stack(
    [np.array(adata.obs['x'].tolist()), np.array(adata.obs['y'].tolist())],
    0
).T

adata.obs_names = "CELL_" + adata.obs_names.astype(str)  # as if otherwise it's a problem for nico
assert np.allclose(
    adata.X.data,
    np.floor(adata.X.data)
)

adata.layers['counts'] = adata.X.copy()

adata_original = adata.copy()


sc.tl.pca(adata)
# sc.pp.neighbors(adata)
sc.pp.neighbors(adata, use_rep='X')  # TODO:provide the tissue section ID
adata.raw = adata

adata_unnorm = adata.copy()


adata_int = sc.read_h5ad(
    args.fname_adata_int
)



# discard the homogenous regions if needed
if args.flag_drop_homogregions:
    adata = adata[
        (adata.obs['region_label'] != "Region5") & (adata.obs['region_label'] != "Region6") 
    ]

if args.flag_drop_homogregions:
    adata_int = adata_int[
        (adata_int.obs['region_label'] != "Region5") & (adata_int.obs['region_label'] != "Region6") 
    ]

if args.flag_drop_homogregions:
    adata_original = adata_original[
        (adata_original.obs['region_label'] != "Region5") & (adata_original.obs['region_label'] != "Region6") 
    ]

if args.flag_drop_homogregions:
    adata_unnorm = adata_unnorm[
        (adata_unnorm.obs['region_label'] != "Region5") & (adata_unnorm.obs['region_label'] != "Region6") 
    ]



# dump the anndata objects, to be read and used by NiCo ===
path_myworkspace = "./NonGit/RunOutputs/{}/".format(str_runname)
os.makedirs(path_myworkspace, exist_ok=True)
path_workspace = path_myworkspace
path_nicoworkspace = path_myworkspace 

adata.write_h5ad(
    os.path.join(
        path_myworkspace,
        'adata.h5ad'
    )
)

adata.write_h5ad(
    os.path.join(
        path_nicoworkspace,
        'adata.h5ad'
    )
)

# sct_spatial

adata.write_h5ad(
    os.path.join(
        path_myworkspace,
        'sct_spatial.h5ad'
    )
)

adata.write_h5ad(
    os.path.join(
        path_nicoworkspace,
        'sct_spatial.h5ad'
    )
)

print("Done!")





annotation_slot = "cell_type" #content_config_ds_annotations_yml[args.ds_name]['obskey_celltype']
showit=True
saveas='png'
transparent_mode = False


# In[ ]:


# parameters of the nico 
inputRadius=0  # ak: recall: Radius=0 means jaxtacrine signalling, 
# annotation_slot='rctd_first_type' #spatial cell type slot, ak: it's the obskey of cell types.
# fname_anndata = ''
do_not_use_following_CT_in_niche=[]



niche_pred_output=sint.spatial_neighborhood_analysis(
    Radius=inputRadius,
    output_nico_dir=path_workspace,
    anndata_object_name='adata.h5ad',
    spatial_cluster_tag=annotation_slot,
    removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

t_begin = time.time()

cov_out=scov.gene_covariation_analysis(
    quepath=path_workspace,
    refpath=path_workspace,
    iNMFmode=False,
    Radius=inputRadius,
    no_of_factors=args.num_factors,
    spatial_integration_modality='single',
    anndata_object_name='adata.h5ad',
    output_niche_prediction_dir=path_workspace,
    ref_cluster_tag=annotation_slot,
    LRdbFilename='NiCoLRdb.txt'
)  # according to gAI, this is related to the NMF step to infer latent factors.

print("Took {} seconds.".format(
    time.time() - t_begin
))


data1=np.load(cov_out.outputname,allow_pickle=True)
data1=data1['weighted_neighborhood_of_factors_in_niche']
data=np.nan_to_num(data1)

featureVector=range(cov_out.no_of_pc,data.shape[1]) # #just neighborhood
X_latent_factors = data[:,0:cov_out.no_of_pc]  # of shape [num_cells x num_factors]
print(X_latent_factors.shape)


dict_ct_to_idxct = {v:k for k, v in niche_pred_output.nameOfCellType.items()}
dict_ct_to_idxct


list_adata_todump = []
for ct in set(adata.obs[annotation_slot].tolist()):
    idx_ct = dict_ct_to_idxct[ct]
    stdnorm_ct, W_ct = cov_out.ADDED_dict_clidsp_to_stdnorm_and_W[str(idx_ct)]  # [num_genes], [num_genes x num_factors]

    # subset the mic part of W_ct
    percent_variance_explained = np.array(cov_out.save_reg_coef[idx_ct][8])  # to be used to split the factors to Xint and Xmic
    W_ct = W_ct[
        :,
        (percent_variance_explained > np.median(percent_variance_explained)).tolist() 
    ]  # [num_genes x num_somefactors]
    if len(W_ct.shape) == 1:
        W_ct = np.expand_dims(W_ct, -1)


    X_latent_factors_ct = X_latent_factors[
        (adata.obs[annotation_slot] == ct).tolist(),
        :
    ]  # [num_cells_CT x num_factors]
    X_latent_factors_ct = X_latent_factors_ct[
        :,
        (percent_variance_explained > np.median(percent_variance_explained)).tolist()
    ]  # [num_cells_CT x num_some_factors]

    if len(X_latent_factors_ct.shape) == 1:
        X_latent_factors_ct = np.expand_dims(X_latent_factors_ct, -1)

    X_mic_ct = X_latent_factors_ct @ (W_ct.T)  # [num_cells_CT x num_genes]

    # un-normalise X_mic_ct
    X_mic_ct = X_mic_ct * np.expand_dims(stdnorm_ct, 0)
    X_mic_ct = np.expm1(X_mic_ct)

    # create adata_ct
    adata_ct = adata_original[
        (adata.obs[annotation_slot] == ct).tolist(),
        :
    ].copy()
    orig_rowsum = np.expand_dims(
        np.array(adata_ct.X.sum(1).tolist()).flatten(),
        -1
    )
    X_mic_ct = X_mic_ct * (orig_rowsum / 10000.0)


    assert isinstance(X_mic_ct, np.ndarray)
    X_mic_ct = X_mic_ct * ((adata_ct.X.toarray() > 0.0) + 0.0)
    X_mic_ct = sparse.csr_matrix(X_mic_ct)
    
    X_mic_ct.data = np.clip(
        X_mic_ct.data,
        a_min=np.zeros_like(adata_ct.X.data),
        a_max=adata_ct.X.data
    )
    X_mic_ct.eliminate_zeros()

    adata_ct.obsm['NiCo_Xmic'] = X_mic_ct
    list_adata_todump.append(adata_ct)


adata_todump = anndata.concat(list_adata_todump)
adata_todump.X = adata_todump.layers['counts'].copy()
adata_todump



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


np_xspl_gt = adata_unnorm.X.toarray() - adata_int.X
dict_output_e = evaluator.eval(
    np_xspl_gt=np_xspl_gt,
    np_xspl_pred=adata_todump.obsm['NiCo_Xmic'],
    np_xobs=adata_unnorm.X.toarray(),
    flag_normalize=False  # OLD: note that for simiv this flag is set to True.
)




dict_todump = {
    'dict_result':dict_output_e
}
with open(os.path.join(path_myworkspace, 'vistoken.pkl'), 'wb') as f:
    pickle.dump(dict_todump, f)
print("Dumped the result.")


# convert `dict_output_e` to a dataframe
path_ouptput_files = path_myworkspace
pattern = r"(\w+)\s+\(among readout == ([\d.]+),\s*total=(\d+)\)"
df_result_finegrained = []
for k in dict_output_e.keys():
    match = re.search(pattern, k)
    assert match
    metric, readout, total = match.groups()
    df_result_finegrained.append([
        'NiCo',
        metric, readout, total, dict_output_e[k]
    ])

df_result_finegrained = pd.DataFrame(
    df_result_finegrained,
    columns=[
        'method',
        'metric', 'readout', 'total', 'value'
    ]
)

df_result_finegrained = df_result_finegrained[df_result_finegrained['metric'] != 'PearsonCorrelation']

df_result_finegrained.to_csv(
    os.path.join(
        path_ouptput_files,
        'df_result_finegrained.csv'
    )
)


# save the overal metric
df_result_overall = []
for metric in set(df_result_finegrained['metric'].tolist()):
    np_numbers = df_result_finegrained[
        df_result_finegrained['metric'] == metric
    ]['value'].tolist()
    np_numbers = np.log(np.array(np_numbers) + 1e-10).mean()


    df_result_overall.append([
        args.fname_adata,
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

