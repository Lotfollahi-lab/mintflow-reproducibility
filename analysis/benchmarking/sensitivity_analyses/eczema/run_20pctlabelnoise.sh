#!/bin/bash
# BSUB -G team361
# BSUB -q gpu-lotfollahi
# BSUB -n 8
# BSUB -M 100000
# BSUB -R "select[mem>100000] rusage[mem=100000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J 20pctlabelnoise
# BSUB -o 20pctlabelnoise_%J.out
# BSUB -e 20pctlabelnoise_%J.err

source /etc/profile.d/modules.sh
module load ISG/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/
python /nfs/team361/sb75/mintflow-reproducibility/analysis/benchmarking/sensitivity_analyses/eczema/nb2_train_20pctlabelnoise.py