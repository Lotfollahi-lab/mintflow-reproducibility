

To run NiCo on synthetic data:
- Install NiCo.
- Clone [NiCo repository](https://github.com/ankitbioinfo/nico_tutorial?tab=readme-ov-file).
- In the cloned NiCo repository, replace the file `NiCo/Covariations.py` with the file `Covariations.py` provided here in the current folder.
- Here in this repository/directory, in `run_NiCo.py` replace the variables `PATH_1` and `PATH_2` to where you have cloned the NiCo repository.
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
