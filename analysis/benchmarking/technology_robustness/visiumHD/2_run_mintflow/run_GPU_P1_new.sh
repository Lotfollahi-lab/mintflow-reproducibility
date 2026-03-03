#!/bin/bash
# BSUB -G team361
# BSUB -q gpu-basement
# BSUB -n 8
# BSUB -M 300000
# BSUB -R "select[mem>300000] rusage[mem=300000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J GPU_P1
# BSUB -o GPU_P1_%J.out
# BSUB -e GPU_P1_%J.err

source ~/.bashrc
module load cellgen/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec11th_mintflowandrapids/
python /nfs/users/nfs_m/ms83/mintflow-revision/nb1_train_P1_7Kgenes.py
