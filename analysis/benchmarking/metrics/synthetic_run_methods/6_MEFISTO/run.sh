

module load ISG/conda

conda activate /nfs/team361/aa36/PythonEnvs_2/env_mefistoV7/

# run on simulated tissue section 1
python run_MEFISTO.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'RunSim1' \
--num_factors_divby2 5 \
--n_inducing 100 \


# run on simulated tissue section 2
python run_MEFISTO.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim2sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim2sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'RunSim2' \
--num_factors_divby2 5 \
--n_inducing 100 \


# run on simulated tissue section 3
python run_MEFISTO.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim3sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim3sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'RunSim3' \
--num_factors_divby2 5 \
--n_inducing 100 \


# run on simulated tissue section 4
python run_MEFISTO.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim4sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim4sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'RunSim4' \
--num_factors_divby2 5 \
--n_inducing 100 \


# run on simulated tissue section 5
python run_MEFISTO.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim5sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim5sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'RunSim5' \
--num_factors_divby2 5 \
--n_inducing 100 \



