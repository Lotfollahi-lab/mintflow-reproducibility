

#!/usr/bin/env python
# coding: utf-8

# This notebook calculates R2 score based on count data probabilities. 

# In[ ]:


import argparse
import types
import os, sys
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
from tqdm.autonotebook import tqdm
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
import pickle

from scvi.distributions import Poisson, NegativeBinomial, ZeroInflatedNegativeBinomial


# In[ ]:


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


# In[ ]:


def row_normalise_ifneeded(X1):

    assert sparse.issparse(X1)
    row_sums = np.array(X1.sum(axis=1)).flatten()

    if (np.max(row_sums) - np.min(row_sums)) >= 1.0:
        return (diags(1.0/row_sums) @ X1) * 10000
    else:
        return X1


# In[ ]:


#parse arguments ========================================
args = types.SimpleNamespace()
if(flag_isrunning_jupyterNB == False):
    parser = argparse.ArgumentParser(description='dsc')


# path_results
if(flag_isrunning_jupyterNB):
    args.path_results = "./NonGit/Slide-seq - kidney/cellcharter/"
else:
    parser.add_argument('--path_results', type=str,
                        help = 'ddd.')

# batch_size
if(flag_isrunning_jupyterNB):
    args.batch_size = 100
else:
    parser.add_argument('--batch_size', type=int,
                        help = 'ddd.')


# num_epochs
if(flag_isrunning_jupyterNB):
    args.num_epochs = 100
else:
    parser.add_argument('--num_epochs', type=int,
                        help = 'ddd.')


if(flag_isrunning_jupyterNB == False):
    args = parser.parse_args()
print("args = {}".format(args)) #======================================================


# In[ ]:


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[ ]:





# # 1. Load gex and embeddings

# In[ ]:


X_gex = sparse.load_npz(
    os.path.join(
        args.path_results,
        'X.npz'
    )
)
X_gex = row_normalise_ifneeded(X_gex)
print(X_gex.shape)


# In[ ]:


for fname_embeddings in os.listdir(args.path_results):
    if fname_embeddings.endswith(".npz"):
        if len(fname_embeddings) >= len("obsm_"):
            if fname_embeddings[0:len('obsm_')] == 'obsm_':
                break  # already found the fname of embeddings

X_embeddings = np.load(
    os.path.join(
        args.path_results,
        fname_embeddings
    )
)['arr_0']
print(X_embeddings.shape)


# In[ ]:





# # 2. Fit the estimator

# In[ ]:


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

        return history_loss, history_perepoch_loss



# In[ ]:


obj_estimator = CalcReconstructionLoglik(
    count_model='Poisson',
    kwargs_decoder={
        'dim_embeddings':X_embeddings.shape[1],
        'dim_hidden':50,
        'num_genes':X_gex.shape[1]
    }
)
obj_estimator.to(device)


# In[ ]:


history_loss, history_perepoch_loss = obj_estimator.fit(
    X_emb=torch.tensor(X_embeddings, device=device).float(),
    X_gex=torch.tensor(X_gex.toarray(), device=device).float(),
    batch_size=args.batch_size,
    num_epochs=args.num_epochs
)


# In[ ]:


history_loss


# In[ ]:


if flag_isrunning_jupyterNB:
    plt.figure()
    plt.plot(
        range(len(history_loss)),
        history_loss
    )
    plt.show()


# In[ ]:


if flag_isrunning_jupyterNB:
    plt.figure()
    plt.plot(
        range(len(history_perepoch_loss)),
        [np.mean(u) for u in history_perepoch_loss]
    )
    plt.show()


# In[ ]:





# # 3. Dump the result

# In[ ]:


path_dump = os.path.join(
    args.path_results,
    'result_overall_NLL/'
)
os.makedirs(path_dump, exist_ok=True)


# In[ ]:


with open(
    os.path.join(
        path_dump,
        'output.pkl'
    ),
    'wb'
) as f:
    pickle.dump(
        [history_perepoch_loss, history_loss],
        f
    )

print("Done!")


# In[ ]:
