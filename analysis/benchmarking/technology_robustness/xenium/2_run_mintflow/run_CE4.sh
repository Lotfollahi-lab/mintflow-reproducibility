#!/bin/bash
# BSUB -G team361
# BSUB -q training-parallel
# BSUB -n 8
# BSUB -M 100000
# BSUB -R "select[mem>100000] rusage[mem=100000]"
# BSUB -gpu "num=1:mode=exclusive_process:gmem=50000"
# BSUB -J CE4
# BSUB -o CE4_%J.out
# BSUB -e CE4_%J.err
# BSUB -U mlteam_burnin_beta

source ~/.bashrc
module load cellgen/conda
conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec11th_mintflowandrapids/
python /nfs/users/nfs_m/ms83/mintflow-revision/psoriasis/nb1_train_CE4-SKI-27-FO-4-S22-A2.py