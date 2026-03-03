#!/bin/bash
# BSUB -G team361
# BSUB -q gpu-lotfollahi
# BSUB -n 8
# BSUB -M 100000
# BSUB -R "select[mem>100000] rusage[mem=100000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=20000"
# BSUB -J Tcells_superficial
# BSUB -o Tcells_superficial_%J.out
# BSUB -e Tcells_superficial_%J.err

source ~/.bashrc
module load cellgen/conda
conda activate /software/cellgen/team361/ms83/envs/mintflow-revision/
python /nfs/users/nfs_m/ms83/mintflow-revision/benchmarking/cellcharter/analysis/R3C1-R3C7_beacon_Tcells_superficial.py
