

module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs_2/env_NiCo/

cd ../../

python run.py \
--fname_anndata /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/4_TimingAndMemoryUsage/0_Data/NonGit/2_ProcessedData/NumGenes_1000/adata_numcells_10000.h5ad \
--num_factors 5 \
--annotation_slot level_2_cell_type \
--path_workspace ./NonGit/Wspace_10000/ \



