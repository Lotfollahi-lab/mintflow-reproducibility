#!/bin/bash
# BSUB -G team361
# BSUB -q training-parallel
# BSUB -n 8
# BSUB -M 200000
# BSUB -R "select[mem>200000] rusage[mem=200000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J lung12
# BSUB -o lung12_%J.out
# BSUB -e lung12_%J.err

source ~/.bashrc
module load cellgen/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/
python /nfs/users/nfs_m/ms83/mintflow-revision/cosmx/nb1_train_lung12.py
