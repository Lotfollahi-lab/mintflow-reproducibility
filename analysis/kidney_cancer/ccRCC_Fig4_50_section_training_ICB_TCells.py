import os
import sys
import gc
import pickle
import time
import glob

import yaml
import mintflow
import mintflow.interface.base_interface

import scanpy as sc
import pandas as pd
import numpy as np
import torch
from tqdm import tqdm
from scipy.stats import wasserstein_distance

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# ── Paths ──
src_dir = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_raw_per_section"
dst_dir = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_raw_per_section_ICB_TC_relabel"
cache_path = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/cache/mintflow_context_tcell_replacement.pkl"
path_output_files = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/ICB_TCell_Training/outputs"
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


# ── ICB T Cell relabelling (skipped if cache already exists) ──
if not os.path.exists(cache_path):
    os.makedirs(dst_dir, exist_ok=True)

    relabel_map = {
        "CD4+ T lymphocyte": "CD4+ T lymphocyte ICB",
        "CD8+ T lymphocyte": "CD8+ T lymphocyte ICB",
        "DN T lymphocyte": "DN T lymphocyte ICB",
    }

    files_to_relabel = sorted([f for f in os.listdir(src_dir) if f.endswith(".h5ad")])
    print(f"Total files to relabel: {len(files_to_relabel)}")

    for i, fname in enumerate(files_to_relabel):
        adata = sc.read_h5ad(os.path.join(src_dir, fname))
        drug_val = adata.obs["drug"].iloc[0]

        if drug_val == "Ipi/Novo":
            new_cats = [v for v in relabel_map.values() if v not in adata.obs["level_3_cell_type"].cat.categories]
            if new_cats:
                adata.obs["level_3_cell_type"] = adata.obs["level_3_cell_type"].cat.add_categories(new_cats)
            mask = adata.obs["level_3_cell_type"].isin(relabel_map.keys())
            adata.obs.loc[mask, "level_3_cell_type"] = adata.obs.loc[mask, "level_3_cell_type"].map(relabel_map)
            adata.obs["level_3_cell_type"] = adata.obs["level_3_cell_type"].cat.remove_unused_categories()
            n_relabeled = mask.sum()
            print(f"[{i+1}/{len(files_to_relabel)}] {fname} - drug=Ipi/Novo, relabeled {n_relabeled} cells")
        else:
            print(f"[{i+1}/{len(files_to_relabel)}] {fname} - drug={drug_val}, no relabeling")

        adata.write_h5ad(os.path.join(dst_dir, fname))
        del adata
        gc.collect()

    print("Relabelling done!")
else:
    print(f"Cache found at {cache_path}; skipping relabelling step.")

# ── Load / build MintFlow context ──
if os.path.exists(cache_path):
    print(f"Loading cached MintFlow context from: {cache_path}")
    with open(cache_path, "rb") as f:
        data_mintflow, dict_all4_configs = pickle.load(f)
else:
    print(f"Building MintFlow context (cache not found at {cache_path})")
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    t0 = time.time()
    data_mintflow, dict_all4_configs = build_mintflow_context(
        dst_dir, num_training_epochs=NUM_TRAINING_EPOCHS
    )
    with open(cache_path, "wb") as f:
        pickle.dump((data_mintflow, dict_all4_configs), f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved context cache: {cache_path} (elapsed {time.time() - t0:.1f}s)")

config_data_train = dict_all4_configs['config_data_train']
config_training = dict_all4_configs['config_training']

gc.collect()

model = mintflow.setup_model(dict_all4_configs=dict_all4_configs, data_mintflow=data_mintflow)
trainer = mintflow.Trainer(dict_all4_configs=dict_all4_configs, model=model, data_mintflow=data_mintflow)

# ── Training loop with evaluation ──
num_epochs = config_training['num_training_epochs']

for index_epoch in tqdm(range(num_epochs), desc='Training epoch'):
    trainer.train_one_epoch()

    # Save model checkpoint
    mintflow.interface.base_interface.dump_model(
        model,
        os.path.join(path_output_files, f"model_epoch_{index_epoch}.pt"),
    )

    # Predict and save for each tissue section
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

    # Evaluate signalling genes for this epoch
    with torch.no_grad():
        df_evaluation_result = mintflow.evaluate_by_known_signalling_genes(
            device=device,
            dict_all4_configs=dict_all4_configs,
            data_mintflow=data_mintflow,
            model=model,
            evalulate_on_sections='all',
            optional_list_colvaltype_toadd=[['training_epoch', index_epoch, 'category']]
        )

    df_evaluation_result.to_pickle(
        os.path.join(path_output_files, f'df_evaluation_result_epoch_{index_epoch}.pkl')
    )

    del df_evaluation_result
    torch.cuda.empty_cache()
    gc.collect()

# ── Post-training analysis ──
pkl_files = sorted(glob.glob(os.path.join(path_output_files, "df_evaluation_result_epoch_*.pkl")))
dfs = []
for f in pkl_files:
    df = pd.read_pickle(f)
    dfs.append(df)

df_toinspect = pd.concat(dfs, ignore_index=True)

df_filtered = df_toinspect[df_toinspect['read_count'] > 30.0]

wd_records = []
for e in range(num_epochs):
    df_epoch = df_filtered[df_filtered['training_epoch'] == e]
    sig = df_epoch.loc[df_epoch['is_among_signalling_genes'] == 'True', 'fraction_assigned_to_Xmic'].values
    nonsig = df_epoch.loc[df_epoch['is_among_signalling_genes'] == 'False', 'fraction_assigned_to_Xmic'].values
    wd = wasserstein_distance(sig, nonsig)
    wd_records.append({'epoch': e, 'wasserstein_distance': wd})
    print(f"Epoch {e}: Wasserstein distance = {wd:.4f}")
df_wd = pd.DataFrame(wd_records)
df_wd.to_csv(os.path.join(path_output_files, 'wasserstein_distances.csv'), index=False)

fig, axes = plt.subplots(1, 2, figsize=(16, 5))
sns.violinplot(
    data=df_filtered, x='training_epoch', y="fraction_assigned_to_Xmic",
    hue="is_among_signalling_genes", ax=axes[0], cut=0
)
axes[0].set_title("ICB_TCell: Disentanglement per epoch")
axes[1].plot(df_wd['epoch'], df_wd['wasserstein_distance'], marker='o')
axes[1].set_xlabel('Epoch')
axes[1].set_ylabel('Wasserstein distance')
axes[1].set_title("ICB_TCell: Wasserstein distance (signalling vs non-signalling)")
axes[1].set_xticks(range(num_epochs))
plt.tight_layout()
fig.savefig(os.path.join(path_output_files, 'analysis_plots.png'), dpi=150)
plt.close(fig)

# Cleanup
del data_mintflow, model, trainer, dfs, df_toinspect, df_filtered
gc.collect()
print("\nCompleted training for ICB_TCell")
