#!/bin/bash
# BSUB -G team361
# BSUB -q training-parallel
# BSUB -n 8
# BSUB -M 100000
# BSUB -R "select[mem>100000] rusage[mem=100000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J 1000genes
# BSUB -o 1000genes_%J.out
# BSUB -e 1000genes_%J.err

source ~/.bashrc
module load ISG/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/
python /nfs/users/nfs_m/ms83/mintflow-revision/old_kidney/nb2_train_1000genes.py