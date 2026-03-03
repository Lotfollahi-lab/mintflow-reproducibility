#!/bin/bash
# BSUB -G team361
# BSUB -q training-parallel
# BSUB -n 8
# BSUB -M 100000
# BSUB -R "select[mem>100000] rusage[mem=100000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=20000"
# BSUB -J TLScore1_500hvg
# BSUB -o TLScore1_500hvg_%J.out
# BSUB -e TLScore1_500hvg_%J.err

source ~/.bashrc
module load cellgen/conda
conda activate /software/cellgen/team361/ms83/envs/mintflow-revision/
python /nfs/users/nfs_m/ms83/mintflow-revision/old_kidney/benchmarking/R1C0_500hvg_TLScore1.py
