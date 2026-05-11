

module load ISG/conda

conda activate /nfs/team361/aa36/PythonEnvs_2/env_simvi/

# run on simulated tissue section 1
python run_SIMVI.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'Sim1' \
--obskey_celltype 'cell_type' \
--neighgraph_num_neighbours 5 \
--num_factors 20 \
--max_epochs 100 \
--batch_size 500 \
--mae_epochs 25 \


# run on simulated tissue section 2
python run_SIMVI.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim2sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim2sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'Sim1' \
--obskey_celltype 'cell_type' \
--neighgraph_num_neighbours 5 \
--num_factors 20 \
--max_epochs 100 \
--batch_size 500 \
--mae_epochs 25 \


# run on simulated tissue section 3
python run_SIMVI.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim3sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim3sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'Sim1' \
--obskey_celltype 'cell_type' \
--neighgraph_num_neighbours 5 \
--num_factors 20 \
--max_epochs 100 \
--batch_size 500 \
--mae_epochs 25 \



# run on simulated tissue section 4
python run_SIMVI.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim4sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim4sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'Sim1' \
--obskey_celltype 'cell_type' \
--neighgraph_num_neighbours 5 \
--num_factors 20 \
--max_epochs 100 \
--batch_size 500 \
--mae_epochs 25 \


# run on simulated tissue section 5
python run_SIMVI.py \
--fname_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim5sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" \
--fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim5sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" \
--flag_drop_homogregions 'False' \
--str_prefix_runname 'Sim1' \
--obskey_celltype 'cell_type' \
--neighgraph_num_neighbours 5 \
--num_factors 20 \
--max_epochs 100 \
--batch_size 500 \
--mae_epochs 25 \


