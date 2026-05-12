


import os, sys
sys.path.append("/nfs/team361/aa36/OnGit/nico_tutorial/NiCo/")
sys.path.append("/nfs/team361/aa36/OnGit/nico_tutorial/")


# if you installed the nico package 


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



# settings ===
# path_workspace = './NonGit/AkOutput_Test_3/'
# fname_anndata = './nico_cerebellum/cerebellum.h5ad'
# annotation_slot='rctd_first_type' #spatial cell type slot, ak: it's the obskey of cell types.

# args ========
parser = argparse.ArgumentParser(
    prog="DataProcessor",
    description="A simple script to process files.",
    epilog="Use -h for more information."
)
parser.add_argument("--path_workspace", type=str, required=True, help="ddd")
parser.add_argument("--fname_anndata", type=str, required=True, help="ddd")
parser.add_argument("--annotation_slot", type=str, required=True, help="ddd")
parser.add_argument("--num_factors", type=int, required=True, help="ddd")
args = parser.parse_args()


path_workspace = args.path_workspace
fname_anndata = args.fname_anndata
annotation_slot = args.annotation_slot  # args =======



# copy the anndata object into the workspace
os.makedirs(
    path_workspace,
    exist_ok=True
)
if not os.path.isfile(
    os.path.join(
        path_workspace,
        fname_anndata.split("/")[-1]
    )
):
    os.system(
        "cp '{}' '{}' ".format(
            fname_anndata,
            path_workspace
        )
    )
    print("Done!")



# parameters of the nico 
inputRadius=0  # ak: recall: Radius=0 means jaxtacrine signalling, 
# annotation_slot='rctd_first_type' #spatial cell type slot, ak: it's the obskey of cell types.
# fname_anndata = ''
do_not_use_following_CT_in_niche=[]

niche_pred_output=sint.spatial_neighborhood_analysis(
    Radius=inputRadius,
    output_nico_dir=path_workspace,
    anndata_object_name=fname_anndata.split("/")[-1],
    spatial_cluster_tag=annotation_slot,
    removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

t_begin = time.time()

cov_out=scov.gene_covariation_analysis(
    Radius=inputRadius,
    no_of_factors=args.num_factors,
    spatial_integration_modality='single',
    anndata_object_name=fname_anndata.split("/")[-1],
    output_niche_prediction_dir=path_workspace,
    ref_cluster_tag=annotation_slot,
    LRdbFilename='NiCoLRdb.txt'
)

print("Took {} seconds.".format(
    time.time() - t_begin
))



