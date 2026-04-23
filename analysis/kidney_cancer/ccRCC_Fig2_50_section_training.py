import os
import sys
import gc
import pickle
import time

import yaml
import mintflow
import mintflow.interface.base_interface

import torch
from tqdm import tqdm

print(sys.executable)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

path_anndata = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_raw_per_section"
cache_path = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/cache/mintflow_context.pkl"
path_output_files = "/nfs/team361/am74/mintflow/revision/train_model_all_data/window_sz_2280_alldatawithbatcheffect/8epochs_save_h5ad/new_hyperparametes/run2/outputs"
os.makedirs(path_output_files, exist_ok=True)

NUM_TRAINING_EPOCHS = 8


def _mintflow_bool_to_str(x):
    if isinstance(x, dict):
        for k, v in x.items():
            if isinstance(v, bool):
                x[k] = "True" if v else "False"
            else:
                _mintflow_bool_to_str(v)
    elif isinstance(x, list):
        for i, v in enumerate(x):
            if isinstance(v, bool):
                x[i] = "True" if v else "False"
            else:
                _mintflow_bool_to_str(v)


def build_mintflow_context(path_anndata, num_training_epochs=8):
    files = sorted([f for f in os.listdir(path_anndata) if f.endswith(".h5ad")])
    if not files:
        raise FileNotFoundError(f"No .h5ad files found in: {path_anndata}")

    num_sections = len(files)
    config_data_train, config_data_evaluation, config_model, config_training = (
        mintflow.get_default_configurations(
            num_tissue_sections_training=num_sections,
            num_tissue_sections_evaluation=num_sections,
        )
    )

    for idx, filename in enumerate(files, start=1):
        key = f"anndata{idx}"
        filepath = os.path.join(path_anndata, filename)

        config_data_train["list_tissue"][key]["file"] = filepath
        config_data_train["list_tissue"][key]["obskey_cell_type"] = "level_3_cell_type"
        config_data_train["list_tissue"][key]["obskey_sliceid_to_checkUnique"] = "donor"
        config_data_train["list_tissue"][key]["obskey_x"] = "x_centroid"
        config_data_train["list_tissue"][key]["obskey_y"] = "y_centroid"
        config_data_train["list_tissue"][key]["obskey_biological_batch_key"] = "batch"
        config_data_train["list_tissue"][key]["config_dataloader_train"]["width_window"] = 1368
        config_data_train["list_tissue"][key]["config_neighbourhood_graph"] = {
            "n_neighs": 10,
            "set_diag": "False",
            "delaunay": "False",
        }

        config_data_evaluation["list_tissue"][key]["file"] = filepath
        config_data_evaluation["list_tissue"][key]["obskey_cell_type"] = "level_3_cell_type"
        config_data_evaluation["list_tissue"][key]["obskey_sliceid_to_checkUnique"] = "donor"
        config_data_evaluation["list_tissue"][key]["obskey_x"] = "x_centroid"
        config_data_evaluation["list_tissue"][key]["obskey_y"] = "y_centroid"
        config_data_evaluation["list_tissue"][key]["obskey_biological_batch_key"] = "batch"
        config_data_evaluation["list_tissue"][key]["config_dataloader_test"]["width_window"] = 1368
        config_data_evaluation["list_tissue"][key]["config_neighbourhood_graph"] = {
            "n_neighs": 10,
            "set_diag": "False",
            "delaunay": "False",
        }

    config_model["coef_xbarint2notbatchID_loss"] = 1.0
    config_model["coef_xbarspl2notbatchID_loss"] = 1.0
    config_model["coef_flowmatchingloss"] = 0.0
    config_model["dict_qname_to_scaleandunweighted"] = (
        "impanddisentgl_int#0.1#True&impanddisentgl_spl#0.0#True"
        "&varphi_enc_int#0.0#True&varphi_enc_spl#0.0#True"
        "&z#0.0#True&sin#0.0#True&sout#0.0#True"
    )
    config_model["coef_loss_CTpredfromZ"] = 100

    config_training["num_training_epochs"] = num_training_epochs
    config_training["flag_use_GPU"] = "True"
    config_training["flag_enable_wandb"] = "False"
    config_training["annealing_decoder_XintXspl_coef_max"] = 0.01

    config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train)
    config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
    config_model = mintflow.verify_and_postprocess_config_model(
        config_model, num_tissue_sections=len(config_data_train)
    )
    config_training = mintflow.verify_and_postprocess_config_training(config_training)

    dict_all4_configs = {
        "config_data_train": config_data_train,
        "config_data_evaluation": config_data_evaluation,
        "config_model": config_model,
        "config_training": config_training,
    }
    _mintflow_bool_to_str(dict_all4_configs)
    data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
    return data_mintflow, dict_all4_configs


if os.path.exists(cache_path):
    print(f"Loading cached MintFlow context from: {cache_path}")
    with open(cache_path, "rb") as f:
        data_mintflow, dict_all4_configs = pickle.load(f)
else:
    print(f"Building MintFlow context (cache not found at {cache_path})")
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    t0 = time.time()
    data_mintflow, dict_all4_configs = build_mintflow_context(
        path_anndata, num_training_epochs=NUM_TRAINING_EPOCHS
    )
    with open(cache_path, "wb") as f:
        pickle.dump((data_mintflow, dict_all4_configs), f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved context cache: {cache_path} (elapsed {time.time() - t0:.1f}s)")

config_data_train = dict_all4_configs['config_data_train']
config_training = dict_all4_configs['config_training']

gc.collect()

model = mintflow.setup_model(
    dict_all4_configs=dict_all4_configs,
    data_mintflow=data_mintflow,
)
trainer = mintflow.Trainer(
    dict_all4_configs=dict_all4_configs,
    model=model,
    data_mintflow=data_mintflow,
)

for index_epoch in tqdm(range(config_training['num_training_epochs']), desc='Training epoch'):
    '''
    IMPORTANT NOTE: To change the number of epochs, set `config_training['num_training_epochs']` in previous cells of this notebook
    and please refrain from changing the for loop here to, e.g., `for index_epoch in tqdm(range(10), ...)`.
    Because MintFlow's annealing module presumes that the number of epochs equals `config_training['num_training_epochs']`.
    '''

    trainer.train_one_epoch()

    mintflow.interface.base_interface.dump_model(
        model,
        os.path.join(path_output_files, f"model_epoch_{index_epoch}.pt"),
    )

    for idx_tissue in range(len(config_data_train)):
        predictions_tissue = mintflow.predict(
            device=device,
            dict_all4_configs=dict_all4_configs,
            data_mintflow=data_mintflow,
            model=model,
            evalulate_on_sections=[idx_tissue],
        )

        adata_current_tissue = (
            data_mintflow['train_list_tissue_section']
            .list_slice[idx_tissue]
            .adata
        )

        for k in predictions_tissue[f'TissueSection {idx_tissue} (zero-based)']:
            adata_current_tissue.obsm[k] = (
                predictions_tissue[f'TissueSection {idx_tissue} (zero-based)'][k]
            )

        adata_current_tissue.write_h5ad(
            os.path.join(
                path_output_files,
                f'adata_epoch_{index_epoch}_tissuesection_{idx_tissue}.h5ad',
            )
        )

        for k in predictions_tissue[f'TissueSection {idx_tissue} (zero-based)']:
            del adata_current_tissue.obsm[k]

        del predictions_tissue
        gc.collect(); gc.collect(); gc.collect()
        print(f"Done for tissue {idx_tissue}")
