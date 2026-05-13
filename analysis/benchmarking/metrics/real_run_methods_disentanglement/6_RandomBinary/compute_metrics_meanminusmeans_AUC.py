

import os, sys
import importlib
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gc
from tqdm.autonotebook import tqdm
from sklearn.metrics import roc_auc_score, f1_score, precision_score, recall_score
from scipy.special import softmax
from typing import List

import torch
import torch.nn as nn
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

# settings ===
batch_size = 1000

class SimpleClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        self.module = nn.Sequential(
            nn.Linear(1, 10),
            nn.ReLU(),
            nn.Linear(10, 2)
        )

    def forward(self, x):
        return self.module(x)
    
ce_loss = nn.CrossEntropyLoss()


for subdir in os.listdir(
    './NonGit/Runs_TisssuesCombined/'
):
    if os.path.isfile(
        os.path.join(
            './NonGit/Runs_TisssuesCombined/',
            subdir,
            'df_metric_one.pkl'
        )
    ):
        df = pd.read_pickle(
            os.path.join(
                './NonGit/Runs_TisssuesCombined/',
                subdir,
                'df_metric_one.pkl'
            )
        )

        df_todump = []

        # compute mean minus mean
        vals_1 = df.query("is_among_signalling_genes==1.0 and read_count>5.0")['fraction_assigned_to_Xmic'].to_numpy()
        vals_2 = df.query("is_among_signalling_genes==0.0 and read_count>5.0")['fraction_assigned_to_Xmic'].to_numpy()
        val_mean_minus_mean = vals_1.mean() - vals_2.mean()
        
        # compute AUC metric
        curr_df = df.query("read_count>5.0")
        ds = torch.utils.data.TensorDataset(
            torch.tensor(curr_df['fraction_assigned_to_Xmic'].to_numpy()).unsqueeze(-1),
            torch.tensor(
                np.array(
                    [int(u) for u in curr_df['is_among_signalling_genes'].tolist()]
                )
            )
        )
        if len(ds) > 0:
            dl_train = torch.utils.data.DataLoader(ds, batch_size=batch_size, shuffle=True)
            dl_test = torch.utils.data.DataLoader(ds, batch_size=batch_size, shuffle=True)
            model = SimpleClassifier()
            model.to(device)
            optim = torch.optim.Adam(model.parameters(), lr=0.001)
    
            for idx_epoch in tqdm(range(10), desc='       Epoch'):
                # train for 1 epoch
                for _, data in enumerate(tqdm(dl_train, desc='         Looping dl_train')):
                    optim.zero_grad()
                    x, y = data
                    loss = ce_loss(model(x.float().to(device)), y.to(device))
                    loss.backward()
                    optim.step()
    
            # compute the metrcis
            with torch.no_grad():
                list_y, list_pred = [], []  # to be used to predict AUC-ROC
                test_loss = 0.0
                for _, data in enumerate(tqdm(dl_test, desc='         Looping dl_train')):
                    x, y = data
                    netout = model(x.float().to(device))
                    loss = ce_loss(netout, y.to(device))
                    test_loss += loss.sum().detach().cpu().numpy()

                    list_y = list_y + y.tolist()
                    list_pred.append(netout.detach().cpu().numpy())
                

                list_pred = np.concat(list_pred, 0)
                
                val_AUC = roc_auc_score(list_y, list_pred[:,1] - list_pred[:,0])

            gc.collect(); gc.collect(); gc.collect()
        
        else:
            val_AUC = np.nan
        

        df_todump.append([
            val_mean_minus_mean, 
            val_AUC
        ])

        df_todump = pd.DataFrame(
            df_todump,
            columns=[
                'mean-minus-mean',
                'AUC'
            ]
        )

        df_todump.to_csv(
            os.path.join(
                './NonGit/Runs_TisssuesCombined/',
                subdir,
                'df_metrics_Fig1h.csv'
            )
        )

        