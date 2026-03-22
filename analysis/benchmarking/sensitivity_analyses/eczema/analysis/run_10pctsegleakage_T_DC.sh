#!/bin/bash
# BSUB -G team361
# BSUB -q gpu-huge
# BSUB -n 8
# BSUB -M 100000
# BSUB -R "select[mem>100000] rusage[mem=100000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=20000"
# BSUB -J T_DC_10pctsegleakage
# BSUB -o T_DC_10pctsegleakage_%J.out
# BSUB -e T_DC_10pctsegleakage_%J.err

source /etc/profile.d/modules.sh
module load cellgen/conda
conda activate /software/cellgen/team361/ms83/envs/mintflow-revision/
python /nfs/team361/sb75/mintflow-reproducibility/analysis/benchmarking/sensitivity_analyses/eczema/analysis/R2C5_10pctsegleakage_T_DC.py
