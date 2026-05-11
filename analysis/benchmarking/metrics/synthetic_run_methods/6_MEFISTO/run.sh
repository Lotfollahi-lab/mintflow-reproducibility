

module load ISG/conda

conda activate /nfs/team361/aa36/PythonEnvs_2/env_mefistoV7/

# run on simulated tissue section 1
python run_MEFISTO.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'TestRun1' \
--num_factors_divby2 5 \
--n_inducing 100 \

