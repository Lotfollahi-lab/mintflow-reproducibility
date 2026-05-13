

#!/usr/bin/env python
# coding: utf-8

# V4 notes: in V3 the NN wasn't set to .eval() mode. V4 corrects that issue.
# V3 notes: The 2-layer NN contains dropout, to avoid overfitting specially for small populations for each cell type.
# V2 notes: R2 is calculated using a 2-layer NN and Poission likelihood. 
# This notebook takes in (i) a dataset name, (ii) obskey_celltype, and - assuming that Step2 NB is run before - computes per cell type R^2.

# In[1]:


import torch
torch.cuda.is_available()


# In[13]:


import os, sys
import scanpy as sc
import yaml
import argparse


from scipy import sparse
import numpy as np
import pandas as pd

from scipy.sparse import diags
from tqdm.autonotebook import tqdm
import gc

import cupy as cp
from cuml.ensemble import RandomForestRegressor as cuRF
from cuml.decomposition import PCA
from cuml.pipeline import Pipeline
from cuml.preprocessing import StandardScaler
from cuml.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor

from datetime import datetime
import pickle

from scvi.distributions import Poisson, NegativeBinomial, ZeroInflatedNegativeBinomial
import torch
import torch.nn as nn
import torch.nn.functional as F

# parse input args ========
parser = argparse.ArgumentParser(description="DDDD")
parser.add_argument("--ds_name", type=str, required=True, help="D sd fs df sdf")
parser.add_argument("--num_epochs", type=int, required=True, help="D sd fs df sdf")
args = parser.parse_args()  # ========


dict_dsname_to_obskeycelltype = {
    'Xenium_old_kidney':'level_2_cell_type',
    'Xenium_old_beacon': 'broad_celltypes',
    'Slide-seq_kidney': 'scANVI_pred',
    'VisiumHD_CRC': 'DeconvolutionLabel1'
}

dict_dsname_to_fnameannda = {
    'Xenium_old_kidney':'/nfs/team361/ms83/data/kidney/xenium_cellcharter/adata_kidney_all_preprocessed_cellcharter_emb.h5ad',
    'Xenium_old_beacon': '/nfs/team361/ms83/data/skin/xenium_cellcharter/adata_beacon_all_preprocessed_cellcharter_emb.h5ad',
    'Slide-seq_kidney': '/nfs/team361/ms83/data/kidney/slideseq_cellcharter/adata_kidney_slideseq_all_preprocessed_cellcharter_emb.h5ad',
    'VisiumHD_CRC': '/nfs/team361/ms83/data/CRC/visiumHD_cellcharter/adata_sc_all_preprocessed_cellcharter_emb.h5ad'
}



# settings ===
ds_name = args.ds_name  # 'VisiumHD - CRC'
obskey_celltype = dict_dsname_to_obskeycelltype[ds_name]  #  'DeconvolutionLabel1'
batch_size = 100
num_epochs = args.num_epochs  #  10
max_numcells_calcR2 = 10000
n_pca_components = 100


# In[4]:


def row_normalise_ifneeded(X1):

    assert sparse.issparse(X1)
    row_sums = np.array(X1.sum(axis=1)).flatten()

    if (np.max(row_sums) - np.min(row_sums)) >= 1.0:
        return (diags(1.0/row_sums) @ X1) * 10000
    else:
        return X1



device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


class CalcReconstructionLoglik(nn.Module):
    def __init__(self, count_model:str, kwargs_decoder:dict):
        super(CalcReconstructionLoglik, self).__init__()

        # check args
        assert isinstance(count_model, str)
        assert count_model in [
            'Poisson',
            'NegativeBinomial',
            'ZeroInflatedNegativeBinomial'
        ]

        # grab args
        self.count_model = count_model

        # make internals
        self.module_decoder = nn.Sequential(
            nn.Linear(kwargs_decoder['dim_embeddings'], kwargs_decoder['dim_hidden']),
            nn.ReLU(),
            nn.Dropout(),
            nn.Linear(kwargs_decoder['dim_hidden'], kwargs_decoder['num_genes']),
            nn.Softmax()
        )

        if count_model == 'Poisson':
            pass
        else:
            self.theta = torch.nn.Parameter(
                torch.empty(
                    size=[kwargs_decoder['num_genes']]
                ).unsqueeze(0),
                requires_grad=True
            )


    def _get_learnable_params(self):
        if self.count_model == 'Poisson':
            return list(self.module_decoder.parameters())
        else:
            return list(self.module_decoder.parameters()) + [self.theta]


    def fit(self, X_emb, X_gex, batch_size:int, num_epochs:int):
        assert isinstance(X_emb, torch.Tensor)
        assert isinstance(X_gex, torch.Tensor)

        # make ds/dl
        ds = torch.utils.data.TensorDataset(
            X_emb,
            X_gex,
            torch.range(0, X_emb.shape[0]-1).unsqueeze(-1)
        )
        dl_train = torch.utils.data.DataLoader(
            ds,
            batch_size=batch_size,
            shuffle=True
        )
        dl_test = torch.utils.data.DataLoader(
            ds,
            batch_size=batch_size,
            shuffle=False
        )

        # the actual fit
        optimiser = torch.optim.Adam(
            self._get_learnable_params(),
            lr=0.0001
        )
        history_loss = []
        history_perepoch_loss = []
        for idx_epoch in range(num_epochs):
            # one epoch
            self.module_decoder.train()
            for _, data in enumerate(tqdm(dl_train, desc='Epoch {}'.format(idx_epoch))):
                optimiser.zero_grad()

                if self.count_model == 'ZeroInflatedNegativeBinomial':
                    loss = ZeroInflatedNegativeBinomial(
                        mu=10000 * self.module_decoder(data[0]),
                        theta=torch.exp(self.theta),
                        zi_logits=torch.logit(torch.tensor(0.0001)).tolist()

                    ).log_prob(data[1])
                elif self.count_model == 'Poisson':
                    loss = Poisson(
                        rate=10000 * self.module_decoder(data[0]),
                        validate_args=False
                    ).log_prob(data[1])


                loss = (-1.0) * loss.sum(1).mean(0)

                loss.backward()
                optimiser.step()

                history_loss.append(
                    loss.detach().cpu().numpy().tolist()
                )

            # evaluate
            with torch.no_grad():
                self.module_decoder.eval()
                loss_this_epoch = []
                for _, data in enumerate(tqdm(dl_test, desc='  eval of epoch {}'.format(idx_epoch))):
                    if self.count_model == 'ZeroInflatedNegativeBinomial':
                        loss = ZeroInflatedNegativeBinomial(
                            mu=10000 * self.module_decoder(data[0]),
                            theta=torch.exp(self.theta),
                            zi_logits=torch.logit(torch.tensor(0.0001)).tolist()

                        ).log_prob(data[1])
                    elif self.count_model == 'Poisson':
                        loss = Poisson(
                            rate=10000 * self.module_decoder(data[0]),
                            validate_args=False
                        ).log_prob(data[1])

                    loss = (-1.0) * loss.sum(1)

                    loss_this_epoch = loss_this_epoch + loss.detach().cpu().numpy().tolist()

                history_perepoch_loss.append(loss_this_epoch)
                self.module_decoder.train()

        return history_loss, history_perepoch_loss




def compute_R2score(X1, X2, n_pca_components):
    assert isinstance(X1, np.ndarray)
    assert isinstance(X2, np.ndarray)

    # list_toret = []  # In V1 it was per-gene R2-score. Now it's a list containing a a per-epoch number.

    obj_estimator = CalcReconstructionLoglik(
        count_model='Poisson',
        kwargs_decoder={
            'dim_embeddings':X1.shape[1],
            'dim_hidden':50,
            'num_genes':X2.shape[1]
        }
    )
    obj_estimator.to(device)

    _, history_perepoch_loss = obj_estimator.fit(
        X_emb=torch.tensor(X1, device=device).float(),
        X_gex=torch.tensor(X2, device=device).float(),
        batch_size=batch_size,
        num_epochs=num_epochs
    )
    
    return [np.mean(u) for u in history_perepoch_loss]



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


# In[ ]:





# In[8]:


dict_methodname_to_dict_ct_to_listR2score = {}

for method_name in [
    'cellcharter',
    'GraphST',
    'MintFlow',
    'NCEM',
    'STAGATE'
]:

    print("Processing for method: {}".format(method_name))

    # load X, embedding, and df_obs =========
    adata_dot_X = sparse.load_npz(
        os.path.join(
            '../Overall_NLL/NonGit',
            ds_name,
            method_name,
            'X.npz'
        )
    )
    adata_dot_X = row_normalise_ifneeded(adata_dot_X)

    df_obs = pd.read_csv(
        os.path.join(
            '../Overall_NLL/NonGit',
            ds_name,
            method_name,
            'obs.csv'
        )
    )

    for fname_embedding in os.listdir(
        os.path.join(
            '../Overall_NLL/NonGit',
            ds_name,
            method_name
        )
    ):
        if len(fname_embedding) >= len("obsm_"):
            if fname_embedding[0:len("obsm_")] == 'obsm_':
                break  # obsm file is found

    np_embeddings = np.load(
        os.path.join(
            '../Overall_NLL/NonGit',
            ds_name,
            method_name,
            fname_embedding
        )
    )['arr_0']

    # loop over cell types
    dict_ct_to_listR2score = {}
    for ct in set(df_obs[obskey_celltype].tolist()):
        X_CT = adata_dot_X[(df_obs[obskey_celltype] == ct).tolist()]
        np_embeddings_CT = np_embeddings[(df_obs[obskey_celltype] == ct).tolist()]

        # subsample cells if needed
        if X_CT.shape[0] >= max_numcells_calcR2:
            list_idx_rowsel = np.random.permutation(X_CT.shape[0]).tolist()[0:max_numcells_calcR2]
            list_idx_rowsel.sort()

            X_CT = X_CT[list_idx_rowsel]
            np_embeddings_CT = np_embeddings_CT[list_idx_rowsel]

        list_r2_score = compute_R2score(
            np_embeddings_CT,
            X_CT.toarray(),
            n_pca_components=n_pca_components
        )

        dict_ct_to_listR2score[ct] = list_r2_score


    dict_methodname_to_dict_ct_to_listR2score[method_name] = dict_ct_to_listR2score

    # clean up
    del adata_dot_X, df_obs, np_embeddings
    gc.collect(); gc.collect(); gc.collect(); gc.collect() 




# In[ ]:





# # Dump the result

# In[ ]:





# In[12]:


str_nbstarttime = "output_created_at_"+datetime.now().strftime("%d/%m/%Y-%H:%M:%S").replace('/', '-').replace(':','-')
str_nbstarttime


# In[ ]:


path_dump = os.path.join(
    './NonGit/Output_perCT_R2scores/',
    ds_name,
    obskey_celltype
)

os.makedirs(path_dump, exist_ok=True)

path_dump = os.path.join(
    path_dump,
    str_nbstarttime + ".pkl"
)


with open(path_dump, 'wb') as f:
    pickle.dump(
        dict_methodname_to_dict_ct_to_listR2score,
        f
    )

print("Done!")

