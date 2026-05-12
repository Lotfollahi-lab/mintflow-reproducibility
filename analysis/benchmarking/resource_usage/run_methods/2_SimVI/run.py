



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

import simvi
from simvi.model import SimVI
import pickle



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
parser.add_argument("--obskey_celltype", type=str, required=True, help="ddd")
parser.add_argument("--neighgraph_num_neighbours", type=int, required=True, help="ddd")
parser.add_argument("--num_factors", type=int, required=True, help="ddd")  # defaults to 20 in SimVI's tutorial.

# training
parser.add_argument("--max_epochs", type=int, required=True, help="ddd")  # 100 in SimVI's tutorial.
parser.add_argument("--batch_size", type=int, required=True, help="ddd")  # 500 in SimVI's tutorial.
parser.add_argument("--mae_epochs", type=int, required=True, help="ddd")  # 25 in SimVI's tutorial.

# parser.add_argument("--annotation_slot", type=str, required=True, help="ddd")
args = parser.parse_args()   # args =======



# load and setup the anndata object
adata_mtg = sc.read_h5ad(args.fname_anndata)
adata_mtg.obs['akrun_simvi_batchID'] = 'batch_1'
adata_mtg.obs['batch'] = 'batch_1'

SimVI.setup_anndata(
    adata_mtg,
    batch_key='batch',
    labels_key=args.obskey_celltype
)
print("Finished setting up the anndata.")


# create the simvi model
edge_index = SimVI.extract_edge_index(
    adata_mtg,
    n_neighbors=args.neighgraph_num_neighbours,
    batch_key='batch'
)

# from pytorch_lightning.utilities.seed import seed_everything
from lightning.pytorch import seed_everything
seed_everything(0)


model = SimVI(
    adata_mtg,
    n_batch=len(set(adata_mtg.obs['batch'])),
    kl_weight=1,
    kl_gatweight=0.01,
    lam_mi=1000,
    permutation_rate=0.5,
    n_spatial=args.num_factors,
    n_intrinsic=args.num_factors
)
print("Created the simvi model.")


print("Started at: {}".format(
    datetime.now()
))
t_tic = time.time()
train_loss, val_loss = model.train(
    edge_index,
    max_epochs=args.max_epochs,
    batch_size=args.batch_size,
    use_gpu=True,
    mae_epochs=args.mae_epochs
)
print("Took {} seconds.".format(time.time() - t_tic))



# evaluation phase ===
z_simvi = model.get_latent_representation(edge_index,representation_kind='intrinsic', give_mean=True)
s_simvi = model.get_latent_representation(edge_index,representation_kind='interaction', give_mean=True)

# spateffect_simvi = model.get_spatial_effect(edge_index,n_neighbors = 50)
adata_mtg.obsm['simvi_z'] = z_simvi + 0.0
adata_mtg.obsm['simvi_s'] = s_simvi + 0.0

se_list, r2_zlist, r2_slist, r2_zpvlist, r2_spvlist, S = model.get_se(
    edge_index,
    adata=adata_mtg,
    num_arch =7,
    Kfold=1,
    transformation='none',
    cell_type_label=args.obskey_celltype,
    batch_label='batch'
)

# dump the results
dict_todump = {
    'se_list':se_list,
    'r2_zlist':r2_zlist,
    'r2_slist':r2_slist,
    'r2_zpvlist':r2_zpvlist, 
    'r2_spvlist':r2_spvlist, 
    'S':S
}


# put the predictions in anndata object to be dumped
for k, v in dict_todump.items():
    if isinstance(v, list):
        assert len(v) == 7
        dict_todump[k] = np.stack(v, -1)

for k, v in dict_todump.items():
    print(
        k,
        type(v),
        v.shape
    )

for k, v in dict_todump.items():
    adata_mtg.uns[k] = v
    
    # if len(v.shape) == 2:
    #     adata_mtg.obsm[k] = v
    # else:
    #     adata_mtg.uns[k] = v



# dump the anndata object containing predictions
os.makedirs(
    args.path_workspace,
    exist_ok=True
)
adata_mtg.write_h5ad(
    os.path.join(
        args.path_workspace,
        'output.h5ad'
    )
)


# with open(
#     os.path.join(
#         args.path_workspace,
#         'output.pkl'
#     ),
#     'wb'
# ) as f:
#     pickle.dump(
#         dict_todump,
#         f
#     )

print("Script finished succesfully!")

# breakpoint()




