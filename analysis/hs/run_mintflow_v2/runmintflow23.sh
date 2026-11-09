#!/bin/bash
# activate conda environment
source ~/.bashrc 
module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec11th_mintflowandrapids/
python /nfs/users/nfs_l/ls34/inflow_notebooks/run_mintflow_v2/run_mintflow23.py
echo "Finished LOADING"


