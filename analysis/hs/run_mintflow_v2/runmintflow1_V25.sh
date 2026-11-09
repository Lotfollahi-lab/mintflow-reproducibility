#!/bin/bash
# activate conda environment
source ~/.bashrc 
module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/ #/nfs/team361/aa36/PythonEnvs_2/env_Dec11th_mintflowandrapids/
python /nfs/users/nfs_l/ls34/inflow_notebooks/run_mintflow_v2/run_mintflow1_v25.py
echo "Finished LOADING"


