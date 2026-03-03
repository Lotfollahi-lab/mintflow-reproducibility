#!/bin/bash
# BSUB -G team361
# BSUB -q training-parallel
# BSUB -n 8
# BSUB -M 200000
# BSUB -R "select[mem>200000] rusage[mem=200000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J puck_01
# BSUB -o puck_01_%J.out
# BSUB -e puck_01_%J.err

source ~/.bashrc
module load cellgen/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/
python /nfs/users/nfs_m/ms83/mintflow-revision/slide-seq/nb1_train_slideseq_01_50epochs.py
