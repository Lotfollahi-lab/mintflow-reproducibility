

module load ISG/conda  

conda activate /nfs/team361/aa36/PythonEnvs_2/env_mefistoV7/

# CRC ===
python run_mefisto.py  \
--ds_name '4_CRC_VisiumHD_Publis' \
--str_prefix_runname 'MEFISTO_CRC_0' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 0 \

# AD ===
python run_mefisto.py  \
--ds_name 'Z_AD_10sections_MintFlowPreprint' \
--str_prefix_runname 'MEFISTO_AD_0' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 0 \


python run_mefisto.py  \
--ds_name 'Z_AD_10sections_MintFlowPreprint' \
--str_prefix_runname 'MEFISTO_AD_1' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 1 \


python run_mefisto.py  \
--ds_name 'Z_AD_10sections_MintFlowPreprint' \
--str_prefix_runname 'MEFISTO_AD_2' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 2 \


python run_mefisto.py  \
--ds_name 'Z_AD_10sections_MintFlowPreprint' \
--str_prefix_runname 'MEFISTO_AD_3' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 3 \


# RCC ===
python run_mefisto.py  \
--ds_name '1024HVGs_1_RCC' \
--str_prefix_runname 'MEFISTO_RCC_0' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 0 \


python run_mefisto.py  \
--ds_name '1024HVGs_1_RCC' \
--str_prefix_runname 'MEFISTO_RCC_1' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 1 \


python run_mefisto.py  \
--ds_name '1024HVGs_1_RCC' \
--str_prefix_runname 'MEFISTO_RCC_2' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 2 \


python run_mefisto.py  \
--ds_name '1024HVGs_1_RCC' \
--str_prefix_runname 'MEFISTO_RCC_3' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 3 \


# KidneySlideSeq
python run_mefisto.py  \
--ds_name 'SlideSeq_Kidney_puck01' \
--str_prefix_runname 'MEFISTO_SlideSeq_Kidney_puck01_0' \
--num_factors_divby2 5 \
--n_inducing 100 \
--index_tissue_section 0 \



python create_allDFs.py


python compute_metrics_meanminusmeans_AUC.py 

