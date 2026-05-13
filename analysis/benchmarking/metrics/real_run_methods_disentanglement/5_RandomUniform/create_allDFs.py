


"""
Creates all metric-1 dataframes for metric1, so the figures can be made from them immediately.
"""

import os, sys
import importlib
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gc

import mintflow
import mintflow.data.for_evaluation.db_signalling_genes
import mintflow.evaluation.base_evaluation

from tqdm.autonotebook import tqdm


# settings ===
obsm_key = 'Xmic_random_uniform'
fname_adata_in_eachsubdir = 'adata_result.h5ad'


def func_make_df_eval(adata_input):

    # make sure adata_input.X contains raw unnormalised counts
    assert np.allclose(
        adata_input.X.data,
        np.floor(adata_input.X.data)
    )
    
    # Load the LRDB
    f = importlib.resources.open_binary(
        "mintflow.data.for_evaluation.db_signalling_genes",
        "df_LRpairs_Armingoletal.txt"
    )

    df_LRpairs = pd.read_csv(f)
    f.close()

    list_known_LRgenes_inDB = [
        genename
        for colname in ['LigName', 'RecName'] for group in df_LRpairs[colname].tolist() for genename in str(group).split("__")
    ]
    list_known_LRgenes_inDB = set(list_known_LRgenes_inDB)  # ~9K genes (all in DB)

    # check the number of adata's genes found in the DB (return if there are None)
    num_found_in_LRDB = len(set(adata_input.var.index.tolist()).intersection(set(list_known_LRgenes_inDB)))
    print("In the gene panel {} genes were found in the list of known signalling genes.".format(num_found_in_LRDB))
    if num_found_in_LRDB == 0:
        return
    
    
    # create the df
    df_toret = []
    for idx_colsel, colsel in enumerate([
        adata_input.var.index.isin(set(list_known_LRgenes_inDB)),
        ~adata_input.var.index.isin(set(list_known_LRgenes_inDB))
    ]):
        X = adata_input[:, colsel].X.copy().toarray()
        np_mask_readcount_gt_zero = X > 0.0
        X_mic_beforescppnormalizetotal = adata_input.obsm[obsm_key][:, colsel].toarray()

        np_read_count = X[np_mask_readcount_gt_zero] + 0.0
        np_count_Xmic = X_mic_beforescppnormalizetotal[np_mask_readcount_gt_zero] + 0.0
        np_fraction_Xmic = np_count_Xmic / np_read_count

        flag_is_among_signalling_genes = [True, False][idx_colsel]
        df_toret.append(
            pd.DataFrame(
                np.stack([
                        np_read_count, np_count_Xmic, np_fraction_Xmic,
                        np.array(np_read_count.shape[0] * [flag_is_among_signalling_genes])
                    ],
                    -1
                ),  # [N x 3]
                columns=[
                    mintflow.evaluation.base_evaluation.EvalDFColname.readcount.value,
                    mintflow.evaluation.base_evaluation.EvalDFColname.count_Xmic.value,
                    mintflow.evaluation.base_evaluation.EvalDFColname.fraction_Xmic.value,
                    mintflow.evaluation.base_evaluation.EvalDFColname.among_signalling_genes.value
                ]
            )
        )

    df_toret = pd.concat(df_toret)

    return df_toret
    
    


for subdir in tqdm(os.listdir("./NonGit/Runs_TisssuesCombined/"), desc="Goind through subdirs"):
    if True: #  str_subdirprefix in subdir:
        print(subdir)

        if os.path.isfile(
            os.path.join(
                "./NonGit/Runs_TisssuesCombined/",
                subdir,
                fname_adata_in_eachsubdir
            )
        ):
            fname_pkl_todump = os.path.join(
                "./NonGit/Runs_TisssuesCombined/",
                subdir,
                'df_metric_one.pkl'
            )

            if os.path.isfile(fname_pkl_todump):
                print("the file already exists")
            else:

                adata = sc.read_h5ad(
                    os.path.join(
                        "./NonGit/Runs_TisssuesCombined/",
                        subdir,
                        fname_adata_in_eachsubdir
                    )
                )

                df_evaluation_result = func_make_df_eval(adata_input=adata)

                df_evaluation_result : pd.DataFrame
                df_evaluation_result.to_pickle(
                    fname_pkl_todump
                )

                # clean up
                del adata
                del df_evaluation_result
                gc.collect(); gc.collect(); gc.collect()
        
            
        else:
            print("File not found!")

