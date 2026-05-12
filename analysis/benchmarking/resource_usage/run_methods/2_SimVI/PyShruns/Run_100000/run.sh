

module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_simvi/

cd ../../

python run.py \
--fname_anndata /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/4_TimingAndMemoryUsage/0_Data/NonGit/2_ProcessedData/NumGenes_1000/adata_numcells_100000.h5ad \
--obskey_celltype level_2_cell_type \
--neighgraph_num_neighbours 5 \
--num_factors 20 \
--max_epochs 100 \
--batch_size 500 \
--mae_epochs 25 \
--path_workspace ./NonGit/Wspace_100000/ \
