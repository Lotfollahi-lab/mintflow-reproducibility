



In this directory there is a separate subdirectory corresponding to each method.


To prepare data: 
- Refer to our **[Data & Weights Registry (Google Sheet)](https://docs.google.com/spreadsheets/d/1MZ575oXGSmGBPvi-qBCIgjDxeWn2tQmCxVbO9LWlhNs/edit?usp=sharing)**.
- Download the files that correspond to the following files
    - VisiumHD
        - /nfs/team361/ms83/data/CRC/visiumHD_mintflow/raw_data/adata_sc_P1_rawformintflow.h5ad
    - AD
        - /nfs/team361/aa36/OnGit/inflow-reproducibility/Analysis/15_10Samples_BEACON/0_PreProcessing/Outputs_Preprocessing/nb1_preprocessing_V1/PreProcessingDatedAt-13-01-2025-15-02-45/BK23_Non_lesional_Baseline.h5ad
        - /nfs/team361/aa36/OnGit/inflow-reproducibility/Analysis/15_10Samples_BEACON/0_PreProcessing/Outputs_Preprocessing/nb1_preprocessing_V1/PreProcessingDatedAt-13-01-2025-15-02-45/BK23_Lesional_Baseline.h5ad
        - /nfs/team361/aa36/OnGit/inflow-reproducibility/Analysis/15_10Samples_BEACON/0_PreProcessing/Outputs_Preprocessing/nb1_preprocessing_V1/PreProcessingDatedAt-13-01-2025-15-02-45/BK30_Day_14.h5ad 
        - /nfs/team361/aa36/OnGit/inflow-reproducibility/Analysis/15_10Samples_BEACON/0_PreProcessing/Outputs_Preprocessing/nb1_preprocessing_V1/PreProcessingDatedAt-13-01-2025-15-02-45/BK22_Lesional_Baseline.h5ad 
    - RCC
        - /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/1_Disentanglement/3_Actual_Benchmarks/1_MintFlow/NonGit/DatasetsWithHVGs/Ds_1_RCC/CV2-KID-0-FT-1.h5ad  
        - /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/1_Disentanglement/3_Actual_Benchmarks/1_MintFlow/NonGit/DatasetsWithHVGs/Ds_1_RCC/CV2-KID-0-FT-2.h5ad  
        - /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/1_Disentanglement/3_Actual_Benchmarks/1_MintFlow/NonGit/DatasetsWithHVGs/Ds_1_RCC/CV2-KID-0-FO-1.h5ad  
        - /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/1_Disentanglement/3_Actual_Benchmarks/1_MintFlow/NonGit/DatasetsWithHVGs/Ds_1_RCC/CV6-KID-0-FT-4.h5ad 
    - SlideSeq
        - /nfs/team361/aa36/OnGit/mintflow-benchmarking/Metrics/1_Disentanglement/3_Actual_Benchmarks/NonGit/ModifiedAnndataObjects/SlideseqKidney/adata_puck_01_rawformintflow_7K_svg.h5ad
 
- In the yaml file `./experiments.yml`, replace the paths to where you have downloaded the corresponding file.

After you downloaded the data, to run each method please refer to each subdirectory and the instructions therein.

