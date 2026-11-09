# To run the script, you need to modify the arguments passed to the .py file (i.e. the args that start with `--`)
# Here is the description of each argument
#  - fname_DB: your database will be created in the path that you specify by this arg. 
#     You can have a single file for all of your tissue sections. 
#  - str_prefix_dataset: for example the short description that you've added in the template google doc.
#  - absfname_anndata: the absoluate path to the anndta of the tissue section.
#  - str_list_windowsize: the set of window sizes (i.e. the 1st variable in README.md) that the script will try.
#  - str_list_pygdltrain_numneighbours: the set of values for the 2nd variable in RAEDME.md that the script will try.
#  - config_neighbourhood_graph: config of the neighbourhood graph. You can leave it untouched.
#  - obskey_x, obskey_y: the column names in adata.obs that contain spatial X and Y coordinates.
#  - obskey_cell_type, obskey_sliceid_to_checkUnique, obskey_biological_batch_key: The args for mintflow that specify the columns
#     of adata.obs containing cell type labels, unique ID of each slice, and biological/technical batch labels.
#    The last 3 arguments don't have any effect on the output of this script.  


module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec11th_mintflowandrapids/

python main_tune_window_size.py \
--fname_DB '/nfs/team361/ls34/MintFlowDBs/dataset_skin_test_windowsizes.db' \
--str_prefix_dataset 'skin_test' \
--absfname_anndata '/nfs/team361/ls34/inflow/final_adatas/dataset2/Melanoma.h5ad' \
--str_list_windowsize "[123, 456]" \
--str_list_pygdltrain_numneighbours "[1, 3]" \
--config_neighbourhood_graph "{'n_neighs': 10, 'set_diag': 'False', 'delaunay': 'False'}" \
--obskey_x 'x_centroid' \
--obskey_y 'y_centroid' \
--obskey_cell_type 'broad_annotation' \
--obskey_sliceid_to_checkUnique 'sample_id' \
--obskey_biological_batch_key 'sample_id' \


 