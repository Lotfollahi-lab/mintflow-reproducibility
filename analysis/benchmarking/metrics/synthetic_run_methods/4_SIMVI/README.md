

To run SIMVI on synthetic data:
- Install SIMVI, v0.1.2.
- Download the following files from our **[Data & Weights Registry (Google Sheet)](https://docs.google.com/spreadsheets/d/1MZ575oXGSmGBPvi-qBCIgjDxeWn2tQmCxVbO9LWlhNs/edit?
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad" 
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim1sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad" 
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim2sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad"
    - fname_intrinsic_part_of_simulated_data "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim2sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad"
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim3sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad"
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim3sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad"
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim4sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad"
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim4sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad"
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim5sim1_2000genes_50000locs_strongincrements/sim1_2000genes_50000locs_strongincrements.h5ad"
    - "/nfs/team361/aa36/OnGit/nichecompass-reproducibility/analysis/data_simulation/WD/sim5sim1_2000genes_50000locs_strongincrements/intrinsicpartof_sim1_2000genes_50000locs_strongincrements.h5ad"

- Correct the paths in `./run.sh` to where you have downloaded the files.
- In `./run.sh`, modify `conda activate` so your actual conda environment is activated.
- Run `./run.sh`.


The metric values will be created and saved in csv files in `./NonGit/RunOutputs/`.
In particular, the results will be in the csv files `df_result_finegrained.csv` and `df_result_overal.csv`.
