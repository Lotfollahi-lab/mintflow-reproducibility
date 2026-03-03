#!/bin/bash
# BSUB -G team361
# BSUB -q training-parallel
# BSUB -n 8
# BSUB -M 300000
# BSUB -R "select[mem>300000] rusage[mem=300000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J GPU_P5
# BSUB -o GPU_P5_%J.out
# BSUB -e GPU_P5_%J.err

source ~/.bashrc
module load cellgen/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/
python /nfs/users/nfs_m/ms83/mintflow-revision/visium_HD/nb1_train_P5_7Kgenes_CTloss.py
