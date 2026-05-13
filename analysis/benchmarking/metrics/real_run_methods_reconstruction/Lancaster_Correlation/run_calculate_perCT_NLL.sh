

module load ISG/conda

conda activate /nfs/team361/aa36/PythonEnvs_2/env-mintflow-feature-recommended-configs/


python calc_perCT_NLL.py \
--ds_name 'Xenium_old_kidney' \
--num_epochs 200 \


python calc_perCT_NLL.py \
--ds_name 'Xenium_old_beacon' \
--num_epochs 200 \

python calc_perCT_NLL.py \
--ds_name 'Slide-seq_kidney' \
--num_epochs 200 \

python calc_perCT_NLL.py \
--ds_name 'VisiumHD_CRC' \
--num_epochs 200 \

