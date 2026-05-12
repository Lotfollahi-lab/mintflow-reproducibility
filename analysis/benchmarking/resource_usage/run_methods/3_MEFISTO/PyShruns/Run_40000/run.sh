

module load cellgen/conda && conda activate /nfs/team361/aa36/PythonEnvs/envmefistov1/

cd ../../

python run.py \
--fname_anndata /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/4_TimingAndMemoryUsage/0_Data/NonGit/2_ProcessedData/NumGenes_1000/adata_numcells_40000.h5ad \
--num_factors_divby2 5 \
--n_inducing 100 \
--path_workspace ./NonGit/Wspace_40000/ \
