
module load ISG/conda  

conda activate /nfs/team361/aa36/PythonEnvs_2/env_NiCo/ 


python run_NiCo.py \
--fname_adata "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_adata_int "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--num_factor 6 \



