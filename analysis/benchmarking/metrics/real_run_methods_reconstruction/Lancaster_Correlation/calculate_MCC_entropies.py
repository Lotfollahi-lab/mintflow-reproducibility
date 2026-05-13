

import os, sys
import scanpy as sc
import yaml
import torch
import argparse
import anndata

import mintflow
import mintflow.evaluation.mcc_entropy

import matplotlib.pyplot as plt
import seaborn as sns

# part arguments ====
parser = argparse.ArgumentParser(description="DDD")
parser.add_argument('--ds_name', type=str, help = 'ddd.')
args = parser.parse_args()


# settings ===
ds_name = args.ds_name #'Xenium old - beacon'

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


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

dict_dsname_to_obskeyslide = {
    'Xenium_old_kidney':'section',
    'Xenium_old_beacon': 'sample_id',
    'Slide-seq_kidney': 'batch',
    'VisiumHD_CRC': 'sample'
}

obskey_celltype = dict_dsname_to_obskeycelltype[ds_name]

adata = sc.read_h5ad(
    dict_dsname_to_fnameannda[ds_name]
)
adata

adata : anndata.AnnData

if 'sample_id' in adata.obs.columns.tolist():
    adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')

    list_bad_sampleID = [
        'output-XETG00272__0027131__Region_2__20240628__124459',
        'output-XETG00155__0043607__BK30-SKI-27-FO-4-S3__20241010__133127',
        'output-XETG00272__0027175__Region_3__20240628__124459',
        'output-XETG00272__0027131__Region_8__20240628__124459'
    ]
    for bad_sampleID in list_bad_sampleID:
        assert adata.obs['sample_id'].value_counts()[bad_sampleID] == 1

    adata = adata[
        ~adata.obs['sample_id'].isin(list_bad_sampleID)
    ]


dict_ct_to_MCCentropy, df_for_vis = mintflow.evaluation.mcc_entropy.get_MCC_entropy(
    adata=adata,
    kwargs_neighbourhood_graph={
        'spatial_key': 'spatial',
        'library_key': dict_dsname_to_obskeyslide[ds_name],
        'set_diag': False,
        'delaunay': False,
        'n_neighs': 10
    },
    obskey_celltype=obskey_celltype,
    device=device,
    k_calcentropy=1,
)

sns.barplot(
    df_for_vis,
    x='cell_type',
    y='MCC_entropy'
)

plt.xticks(rotation=90)
plt.show()


# Dump the result
path_dump = './NonGit/PrecomputedMCCEntropies/{}/{}/'.format(
    ds_name,
    obskey_celltype
)
os.makedirs(path_dump, exist_ok=True)

df_for_vis.to_csv(
    os.path.join(
        path_dump,
        'df.csv'
    )
)

sns.barplot(
    df_for_vis,
    x='cell_type',
    y='MCC_entropy'
)

plt.xticks(rotation=90)
plt.savefig(
    os.path.join(
        path_dump,
        'log.png'
    ),
    bbox_inches='tight',
    pad_inches=0
)

print("Done!")
