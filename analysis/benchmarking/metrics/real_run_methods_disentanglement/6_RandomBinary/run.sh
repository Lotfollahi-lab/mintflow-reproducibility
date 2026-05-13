


module load ISG/conda  

conda activate /nfs/team361/aa36/PythonEnvs_2/env-mintflow-feature-recommended-configs/ 


# run on dataset 1
python run_RandomBinary.py \
--ds_name '4_CRC_VisiumHD_Publis' \
--str_prefix_runname 'RandomBinary_CRC' \


# run on dataset 2
python run_RandomBinary.py \
--ds_name 'Z_AD_10sections_MintFlowPreprint' \
--str_prefix_runname 'RandomBinary_AD' \

# run on dataset 3
python run_RandomBinary.py \
--ds_name '1024HVGs_1_RCC' \
--str_prefix_runname 'RandomBinary_RCC' \

# run on dataset 4
python run_RandomBinary.py \
--ds_name 'SlideSeq_Kidney_puck01' \
--str_prefix_runname 'RandomBinary_KidneySlideSeq' \



python create_allDFs.py

python compute_metrics_meanminusmeans_AUC.py 


