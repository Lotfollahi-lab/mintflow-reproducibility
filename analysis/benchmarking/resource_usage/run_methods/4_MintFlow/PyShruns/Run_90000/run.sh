

module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_Dec27th_mintflow/

cd ../../

python run.py \
--fname_anndata /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/4_TimingAndMemoryUsage/0_Data/NonGit/2_ProcessedData/NumGenes_1000/adata_numcells_90000.h5ad \
--obskey_cell_type level_2_cell_type \
--dim_sz 100 \
--path_workspace ./NonGit/Wspace_90000/ \
--width_window 900 \
--wandb_runname MintFLow_Timing_90K \
