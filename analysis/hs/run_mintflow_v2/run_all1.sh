#!/bin/bash
module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec11th_mintflowandrapids/

ADATAS_DIR="/nfs/team361/ls34/inflow/final_adatas/dataset1"

for f in "$ADATAS_DIR"/*.h5ad; do
    python main_tune_window_size.py \
        --fname_DB '/nfs/team361/ls34/MintFlowDBs/dataset_skin1all_windowsizes.db' \
        --str_prefix_dataset 'skin1_all' \
        --absfname_anndata "$f" \
        --str_list_windowsize "[123, 456]" \
        --str_list_pygdltrain_numneighbours "[1, 3]" \
        --config_neighbourhood_graph "{'n_neighs': 10, 'set_diag': 'False', 'delaunay': 'False'}" \
        --obskey_x 'x_centroid' \
        --obskey_y 'y_centroid' \
        --obskey_cell_type 'broad_annotation' \
        --obskey_sliceid_to_checkUnique 'sample_id' \
        --obskey_biological_batch_key 'sample_id'
    sleep 500
done


 