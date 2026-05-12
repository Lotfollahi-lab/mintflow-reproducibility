



import numpy as np

if not hasattr(np, 'mat'):
    np.mat = np.asmatrix

import os, sys
import matplotlib.pyplot as plt 
from matplotlib.collections import PatchCollection
import time
import scanpy as sc
import argparse
from datetime import datetime
from mofapy2.run.entry_point import entry_point
import pickle
import scipy

import torch
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)



# args ========
parser = argparse.ArgumentParser(
    prog="DataProcessor",
    description="A simple script to process files.",
    epilog="Use -h for more information."
)

# general
parser.add_argument("--path_workspace", type=str, required=True, help="ddd")
parser.add_argument("--fname_anndata", type=str, required=True, help="ddd")

# model
parser.add_argument("--num_factors_divby2", type=int, required=True, help="ddd")  # divby2, so it beomces 'intrinsic + microenv' 

# training
parser.add_argument("--n_inducing", type=int, required=True, help='ddd')  # DDD

# parser.add_argument("--annotation_slot", type=str, required=True, help="ddd")
args = parser.parse_args()   # args =======



# load and setup the anndata object
adata = sc.read_h5ad(args.fname_anndata)  # adata.X contains the raw counts at this point


sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)


# breakpoint()


# setup t sdsfdfs
ent = entry_point()
ent.set_data_options(use_float32=True)
ent.set_data_from_anndata(adata)#, features_subset="highly_variable")
ent.set_model_options(factors=2*args.num_factors_divby2)


ent.set_train_options(gpu_device=device)

n_inducing = args.n_inducing# 1000

ent.set_covariates([adata.obsm["spatial"]], covariates_names=["x_centroid", "y_centroid"])
ent.set_smooth_options(sparseGP=True, frac_inducing=n_inducing/adata.n_obs,
                       start_opt=10, opt_freq=10)


# as if scipy.empty does not exist anymore, even in very old versions of scipy.
scipy.empty = np.empty



ent.build()
t_begin = time.time()
ent.run()
print("Training took {} seconds.".format(time.time() - t_begin)) # with 22 factors, training took 193 secs.


# dump the results

np_factors = ent.model.nodes["Z"].getExpectations()["E"]
np_weights = ent.model.nodes["W"].getExpectations()[0]["E"]
np_scales_touse = ent.model.train_stats['scales']


dict_todump = {
    'np_factors':np_factors,
    'np_weights':np_weights,
    'np_scales_touse':np_scales_touse
}


os.makedirs(
    args.path_workspace,
    exist_ok=True
)

adata.obsm['MEFISTO_np_factors'] = np_factors
adata.uns['MEFISTO_np_weights'] = np_weights
adata.uns['MEFISTO_np_scales_touse'] = np_scales_touse

adata.write_h5ad(
    os.path.join(
        args.path_workspace,
        'output.h5ad'
    )
)

print("Script finished succesfully!")

# breakpoint()



