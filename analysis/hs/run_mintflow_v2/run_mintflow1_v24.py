#!/usr/bin/env python
# coding: utf-8
# %%

# %%

path_output_files = "/lustre/scratch126/cellgen/lotfollahi/ls34/inflow/outputcheckpoints4/"


import wandb
wandb.login(key='9ac9ae031080bf037dc34bbe2fd3aaa9c89d5c60')
wandb.init(entity='lotfollahi')

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


# %%


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# %%


BASE_DIR =  '/nfs/team361/ls34/inflow/final_adatas/dataset1_v2/'
import os
os.listdir(BASE_DIR)


# %%


len(os.listdir(BASE_DIR))


# %%


path_anndata = '/nfs/team361/ls34/inflow/outputs/dataset1/'


# %%
CELLTYPEKEY= "broad_celltypes"



for i, index_tissue_section in enumerate(['output-XETG00155__0032945__AX34-SKI-0-FO-1-S3-A1__20241217__145711.h5ad',
 'output-XETG00272__0027131__Region_7__20240628__124459.h5ad',
 'output-XETG00272__0021042__AX29-SKI-0-FO-2-S4__20240815__114049.h5ad',
 'output-XETG00272__0027175__Region_2__20240628__124459.h5ad',
 'output-XETG00155__0032957__AX38-SKI-0-FO-1-S3-B1__20241217__145711.h5ad',
 'output-XETG00155__0043607__BK30-SKI-28-FO-1-S3-4__20241010__133127.h5ad',
 'output-XETG00272__0021042__AX30-SKI-0-FO-1-S3__20240815__114049.h5ad',
 'output-XETG00272__0021042__AX30-SKI-0-FO-2-S3__20240815__114049.h5ad',
 'output-XETG00155__0043607__BK30-SKI-27-FO-1-2-S3__20241010__133127.h5ad',
 'output-XETG00272__0027175__Region_7__20240628__124459.h5ad',
 'output-XETG00155__0032957__AX37-SKI-0-FO-2-S3-A2__20241217__145711.h5ad',
 'output-XETG00272__0041416__BK23-SKI-27-FO-1-S6-C2__20241107__113553.h5ad',
 'output-XETG00155__0043607__BK39-SKI-27-FO-2-S3-4__20241010__133127.h5ad',
 'output-XETG00272__0021042__AX29-SKI-0-FO-1-S4__20240815__114049.h5ad',
 'output-XETG00272__0041416__BK22-SKI-27-FO-1-S10-A1__20241107__113553.h5ad',
 'output-XETG00155__0032945__AX34-SKI-0-FO-2-S3-A2__20241217__145711.h5ad',
 'output-XETG00272__0041416__BK23-SKI-27-FO-2-S6-D1_BK23-SKI-27-FO-4-S6-E1__20241107__113553.h5ad',
 'output-XETG00155__0043607__BK39-SKI-27-FO-1-S3__20241010__133127.h5ad']):
    adata = sc.read_h5ad(BASE_DIR + index_tissue_section    )
    adata.obs['TissueSectionID_for_MintFlow'] = index_tissue_section  # the unique ID assigned to this tissue section
    adata.obs['batchID_for_MintFlow'] = index_tissue_section  # i.e. each tissue section is assumed to come from a different batch

    # convert the added column to type "category" (to avoid the slow down caused by having integer columns in pandas dataframes)
    adata.obs['TissueSectionID_for_MintFlow'] = adata.obs['TissueSectionID_for_MintFlow'].astype("category")
    adata.obs['batchID_for_MintFlow'] = adata.obs['batchID_for_MintFlow'].astype("category")

    # save the anndata object back to disk
    adata.write_h5ad(
        f'{path_anndata}/tissue_section_forMintFlow_{i+1}.h5ad'
    )


# %%


config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
    num_tissue_sections_training=len(os.listdir(BASE_DIR)),
    num_tissue_sections_evaluation=len(os.listdir(BASE_DIR))
)


# %%


WINDOW_SIZE=500


# %%


config_data_train['list_tissue']['anndata1']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_1.h5ad')
#   the absolute path to anndata object of tissue section 1 on disk.


config_data_train['list_tissue']['anndata1']['obskey_cell_type'] = CELLTYPEKEY
#   meaning that for the 1st tissue section, cell type labels are provided in `broad_celltypes` column of `adata.obs`.


config_data_train['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
#   meaning that for the 1st tissue section, tissue section ID (i.e. slice ID) is provided in `info_id` column of `adata.obs`


config_data_train['list_tissue']['anndata1']['obskey_x'] = 'x_centroid'
#   meaning that for the 1st tissue section, spatial x coordinates are provided in `x_centroid` column of `adata.obs`


config_data_train['list_tissue']['anndata1']['obskey_y'] = 'y_centroid'
#   meaning that for the 1st tissue section, spatial y coordinates are provided in `y_centroid` column of `adata.obs`



config_data_train['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
#   meaning that for the 1st tissue section, batch identifier is provided in `info_id` column of `adata.obs`


config_data_train['list_tissue']['anndata1']['config_dataloader_train']['width_window'] = WINDOW_SIZE
#   For tissue section one, the widow size of tissue crops used by the customised dataloader 
#      (the dataloader desribed in Supplementary Fig. 16 of paper).
#   The larger this number, the larger the tissue crops, and the bigger the subset of cells in each training iteration.
#      This implies that more GPU memory would be required during training.
#   In this notebook after calling `mintflow.setup_data` in Sec 6 the crop(s) are shown on tissue, 
#      with some information on image titles which can help you tune this parameter.
#   In the manuscript we used `width_window` values between 300 and 800 depending on dataset.


config_data_train['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}
# 


# %%
CELLTYPEKEY='broad_celltypes'
#"fine_annotation"

config_data_train['list_tissue']['anndata2']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_2.h5ad')
config_data_train['list_tissue']['anndata2']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata2']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata2']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata2']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata2']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata2']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata2']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_train['list_tissue']['anndata3']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_3.h5ad')
config_data_train['list_tissue']['anndata3']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata3']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata3']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata3']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata3']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata3']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata3']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_train['list_tissue']['anndata4']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_4.h5ad')
config_data_train['list_tissue']['anndata4']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata4']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata4']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata4']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata4']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata4']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata4']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_train['list_tissue']['anndata5']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_5.h5ad')
config_data_train['list_tissue']['anndata5']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata5']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata5']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata5']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata5']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata5']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata5']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_train['list_tissue']['anndata6']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_6.h5ad')
config_data_train['list_tissue']['anndata6']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata6']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata6']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata6']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata6']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata6']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata6']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_train['list_tissue']['anndata7']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_7.h5ad')
config_data_train['list_tissue']['anndata7']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata7']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata7']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata7']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata7']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata7']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata7']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_train['list_tissue']['anndata8']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_8.h5ad')
config_data_train['list_tissue']['anndata8']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata8']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata8']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata8']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata8']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata8']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata8']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}




config_data_train['list_tissue']['anndata9']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_9.h5ad')
config_data_train['list_tissue']['anndata9']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata9']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata9']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata9']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata9']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata9']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata9']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# ---------------- anndata10 ----------------
config_data_train['list_tissue']['anndata10']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_10.h5ad')
config_data_train['list_tissue']['anndata10']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata10']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata10']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata10']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata10']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata10']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata10']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata11 ----------------
config_data_train['list_tissue']['anndata11']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_11.h5ad')
config_data_train['list_tissue']['anndata11']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata11']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata11']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata11']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata11']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata11']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata11']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata12 ----------------
config_data_train['list_tissue']['anndata12']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_12.h5ad')
config_data_train['list_tissue']['anndata12']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata12']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata12']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata12']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata12']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata12']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata12']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata13 ----------------
config_data_train['list_tissue']['anndata13']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_13.h5ad')
config_data_train['list_tissue']['anndata13']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata13']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata13']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata13']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata13']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata13']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata13']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata14 ----------------
config_data_train['list_tissue']['anndata14']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_14.h5ad')
config_data_train['list_tissue']['anndata14']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata14']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata14']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata14']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata14']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata14']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata14']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata15 ----------------
config_data_train['list_tissue']['anndata15']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_15.h5ad')
config_data_train['list_tissue']['anndata15']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata15']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata15']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata15']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata15']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata15']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata15']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata16 ----------------
config_data_train['list_tissue']['anndata16']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_16.h5ad')
config_data_train['list_tissue']['anndata16']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata16']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata16']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata16']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata16']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata16']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata16']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata17 ----------------
config_data_train['list_tissue']['anndata17']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_17.h5ad')
config_data_train['list_tissue']['anndata17']['obskey_cell_type'] = CELLTYPEKEY
config_data_train['list_tissue']['anndata17']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata17']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata17']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata17']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata17']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata17']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata18 ----------------
config_data_train['list_tissue']['anndata18']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_18.h5ad')
config_data_train['list_tissue']['anndata18']['obskey_cell_type'] = CELLTYPEKEY 
config_data_train['list_tissue']['anndata18']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_train['list_tissue']['anndata18']['obskey_x'] = 'x_centroid'
config_data_train['list_tissue']['anndata18']['obskey_y'] = 'y_centroid'
config_data_train['list_tissue']['anndata18']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_train['list_tissue']['anndata18']['config_dataloader_train']['width_window'] = WINDOW_SIZE
config_data_train['list_tissue']['anndata18']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}
 

    
    
    
    
    
    
    
    
    
    
    
    
    
    







# %%


config_data_evaluation['list_tissue']['anndata1']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_1.h5ad')
config_data_evaluation['list_tissue']['anndata1']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata1']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata1']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata1']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata1']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata1']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata1']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata2']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_2.h5ad')
config_data_evaluation['list_tissue']['anndata2']['obskey_cell_type'] =  CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata2']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata2']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata2']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata2']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata2']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata2']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata3']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_3.h5ad')
config_data_evaluation['list_tissue']['anndata3']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata3']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata3']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata3']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata3']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata3']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata3']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata4']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_4.h5ad')
config_data_evaluation['list_tissue']['anndata4']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata4']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata4']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata4']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata4']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata4']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata4']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata5']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_5.h5ad')
config_data_evaluation['list_tissue']['anndata5']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata5']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata5']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata5']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata5']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata5']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata5']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata6']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_6.h5ad')
config_data_evaluation['list_tissue']['anndata6']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata6']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata6']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata6']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata6']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata6']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata6']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata7']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_7.h5ad')
config_data_evaluation['list_tissue']['anndata7']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata7']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata7']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata7']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata7']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata7']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata7']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


# %%


config_data_evaluation['list_tissue']['anndata8']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_8.h5ad')
config_data_evaluation['list_tissue']['anndata8']['obskey_cell_type'] =CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata8']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata8']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata8']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata8']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata8']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata8']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}


config_data_evaluation['list_tissue']['anndata9']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_9.h5ad')
config_data_evaluation['list_tissue']['anndata9']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata9']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata9']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata9']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata9']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata9']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata9']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata10 ----------------
config_data_evaluation['list_tissue']['anndata10']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_10.h5ad')
config_data_evaluation['list_tissue']['anndata10']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata10']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata10']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata10']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata10']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata10']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata10']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata11 ----------------
config_data_evaluation['list_tissue']['anndata11']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_11.h5ad')
config_data_evaluation['list_tissue']['anndata11']['obskey_cell_type'] =  CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata11']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata11']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata11']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata11']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata11']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata11']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata12 ----------------
config_data_evaluation['list_tissue']['anndata12']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_12.h5ad')
config_data_evaluation['list_tissue']['anndata12']['obskey_cell_type'] =  CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata12']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata12']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata12']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata12']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata12']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata12']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata13 ----------------
config_data_evaluation['list_tissue']['anndata13']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_13.h5ad')
config_data_evaluation['list_tissue']['anndata13']['obskey_cell_type'] =  CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata13']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata13']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata13']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata13']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata13']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata13']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata14 ----------------
config_data_evaluation['list_tissue']['anndata14']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_14.h5ad')
config_data_evaluation['list_tissue']['anndata14']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata14']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata14']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata14']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata14']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata14']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata14']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata15 ----------------
config_data_evaluation['list_tissue']['anndata15']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_15.h5ad')
config_data_evaluation['list_tissue']['anndata15']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata15']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata15']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata15']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata15']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata15']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata15']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata16 ----------------
config_data_evaluation['list_tissue']['anndata16']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_16.h5ad')
config_data_evaluation['list_tissue']['anndata16']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata16']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata16']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata16']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata16']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata16']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata16']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata17 ----------------
config_data_evaluation['list_tissue']['anndata17']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_17.h5ad')
config_data_evaluation['list_tissue']['anndata17']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata17']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata17']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata17']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata17']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata17']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata17']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}

# ---------------- anndata18 ----------------
config_data_evaluation['list_tissue']['anndata18']['file'] = os.path.join(path_anndata, 'tissue_section_forMintFlow_18.h5ad')
config_data_evaluation['list_tissue']['anndata18']['obskey_cell_type'] = CELLTYPEKEY
config_data_evaluation['list_tissue']['anndata18']['obskey_sliceid_to_checkUnique'] = 'TissueSectionID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata18']['obskey_x'] = 'x_centroid'
config_data_evaluation['list_tissue']['anndata18']['obskey_y'] = 'y_centroid'
config_data_evaluation['list_tissue']['anndata18']['obskey_biological_batch_key'] = 'batchID_for_MintFlow'
config_data_evaluation['list_tissue']['anndata18']['config_dataloader_test']['width_window'] = WINDOW_SIZE
config_data_evaluation['list_tissue']['anndata18']['config_neighbourhood_graph'] = {
    'n_neighs': 10,
    'set_diag': 'False',
    'delaunay': 'False',
}













































# %%


config_model['coef_xbarint2notbatchID_loss'] = 1.0
config_model['coef_xbarspl2notbatchID_loss'] = 1.0


# %%


config_model['coef_flowmatchingloss'] = 0.0
config_model['dict_qname_to_scaleandunweighted'] = "impanddisentgl_int#0.1#True&impanddisentgl_spl#0.0#True&varphi_enc_int#0.0#True&varphi_enc_spl#0.0#True&z#0.0#True&sin#0.0#True&sout#0.0#True"


# %%


config_training['num_training_epochs'] = 20
# number of training epochs, i.e. the number of times the model sees the dataset during training.

config_training['flag_use_GPU'] = 'True'
# whether GPU is used.

config_training['flag_enable_wandb'] = 'True'
# if set to True, during training different loss terms are logged to wandb.
# It's highly recommended to enable wandb. Please refer to wandb website for more info: `wandb.ai`



config_training['wandb_project_name'] = 'MintFlow'
# wandb project name (ignored if `config_training['flag_enable_wandb']` is set to False)


# %% changing settings
config_training['wandb_run_name'] = 'Mintflow_inflammation4'
# wandb run name (ignored if `config_training['flag_enable_wandb']` is set to False)
config_training['annealing_decoder_XintXspl_coef_max']=0.01
config_model['coef_loss_CTpredfromZ']=1.0





# %%


config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train) 


# %%


config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)


# %%

config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))  


config_training = mintflow.verify_and_postprocess_config_training(config_training) 


# %%


print("Finished verifying the 4 configuration objects.")


# %%


dict_all4_configs = {
    'config_data_train':config_data_train,
    'config_data_evaluation':config_data_evaluation,
    'config_model':config_model,
    'config_training':config_training
}


# %%


data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)


# %%


model = mintflow.setup_model(
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow
)


# %%


trainer = mintflow.Trainer(
    dict_all4_configs=dict_all4_configs,
    model=model,
    data_mintflow=data_mintflow
)


# %%


# TODO:MODIFY: the path where che


# %%


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
    # df_evaluation_result = mintflow.evaluate_by_known_signalling_genes(
    #     device=device,
    #     dict_all4_configs=dict_all4_configs,
    #     data_mintflow=data_mintflow,
    #     model=model,
    #     evalulate_on_sections='all',
    #     optional_list_colvaltype_toadd=[['training_epoch', index_epoch, 'category']]
    # )
    # df_evaluation_result.to_pickle(
    #     os.path.join(
    #         path_output_files,
    #         'df_evaluation_result_epoch_{}.pkl'.format(index_epoch)
    #     )
    # )

    # save the checkpoint
    # mintflow.dump_checkpoint(
    #     model=model,
    #     data_mintflow=data_mintflow,
    #     dict_all4_configs=dict_all4_configs,
    #     path_dump=os.path.join(path_output_files, "checkpoint_epoch_{}.pt".format(index_epoch)),
    # ) 
    if index_epoch >= (config_training['num_training_epochs'] - 1 - 3): # i.e., if it's the last 3 epochs.
        # save the checkpoint
        mintflow.dump_checkpoint(
            model=model,
            data_mintflow=data_mintflow,
            dict_all4_configs=dict_all4_configs,
            path_dump=os.path.join(path_output_files, "checkpoint_epoch_{}.pt".format(index_epoch)),
        )

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


# %%





# %%





# %%




