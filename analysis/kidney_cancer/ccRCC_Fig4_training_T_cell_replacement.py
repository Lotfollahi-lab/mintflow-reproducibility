import os
import gc
import pickle
import time
import glob

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

run_label = "unseen_scRNAseq_subtype"
dst_dir = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Revisions/adata_raw_per_section_unseen_scRNAseq_subtype"
cache_path = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/cache/mintflow_context_tcell_replacement_unseen_scRNAseq_subtype.pkl"
path_output_files = "/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/ICB_TCell_Training_unseen_scRNAseq_subtype/outputs"
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


def _restore_bools(x):
    # inverse of _mintflow_bool_to_str: verify_and_postprocess converts flag_* strings to
    # real bools, but the cached/stringified configs turn them back to "True"/"False".
    # setup_model / Trainer on this mintflow branch assert real bools, so restore them.
    if isinstance(x, dict):
        for k, v in x.items():
            if v == "True":
                x[k] = True
            elif v == "False":
                x[k] = False
            else:
                _restore_bools(v)
    elif isinstance(x, list):
        for i, v in enumerate(x):
            if v == "True":
                x[i] = True
            elif v == "False":
                x[i] = False
            else:
                _restore_bools(v)


def build_dict_all4_configs(path_anndata, num_training_epochs=8):
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

    return {
        "config_data_train": config_data_train,
        "config_data_evaluation": config_data_evaluation,
        "config_model": config_model,
        "config_training": config_training,
    }


def build_mintflow_context(path_anndata, num_training_epochs=8):
    dict_all4_configs = build_dict_all4_configs(path_anndata, num_training_epochs)
    _mintflow_bool_to_str(dict_all4_configs)
    data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
    return data_mintflow, dict_all4_configs


# Post-verify form with real bools, as setup_model and generate_insilico_ST_data consume it:
# setup_data mutates dict_all4_configs to the stringified form, and the cache branch skips the build.
path_dictall4configs = os.path.join(path_output_files, "dict_all4_configs.pkl")
if not os.path.exists(path_dictall4configs):
    with open(path_dictall4configs, "wb") as f:
        pickle.dump(build_dict_all4_configs(dst_dir, NUM_TRAINING_EPOCHS), f)
    print(f"Wrote {path_dictall4configs}")

# ── Load / build MintFlow context ──
# dst_dir is pre-built by the companion notebook (treatment-naive CV sections with the
# top-5%-per-subtype ICB-scoring T cells relabelled to CD4+/CD8+/DN T lymphocyte ICB).
if os.path.exists(cache_path):
    print(f"Loading cached MintFlow context from: {cache_path}")
    with open(cache_path, "rb") as f:
        data_mintflow, dict_all4_configs = pickle.load(f)
else:
    if not (os.path.isdir(dst_dir) and any(f.endswith(".h5ad") for f in os.listdir(dst_dir))):
        raise FileNotFoundError(
            f"No .h5ad files in {dst_dir}. Run the companion notebook first to build it."
        )
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

# config_model / config_training must carry real bools for setup_model / Trainer
# (config_data stays stringified for setup_data's internal parsing).
_restore_bools(dict_all4_configs['config_model'])
_restore_bools(dict_all4_configs['config_training'])
dict_all4_configs['config_training']['num_training_epochs'] = NUM_TRAINING_EPOCHS

model = mintflow.setup_model(dict_all4_configs=dict_all4_configs, data_mintflow=data_mintflow)
trainer = mintflow.Trainer(dict_all4_configs=dict_all4_configs, model=model, data_mintflow=data_mintflow)

# ── Training loop with evaluation ──
num_epochs = config_training['num_training_epochs']

for index_epoch in tqdm(range(num_epochs), desc='Training epoch'):
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
dfs = [pd.read_pickle(f) for f in pkl_files]
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
axes[0].set_title(f"{run_label}: Disentanglement per epoch")
axes[1].plot(df_wd['epoch'], df_wd['wasserstein_distance'], marker='o')
axes[1].set_xlabel('Epoch')
axes[1].set_ylabel('Wasserstein distance')
axes[1].set_title(f"{run_label}: Wasserstein distance (signalling vs non-signalling)")
axes[1].set_xticks(range(num_epochs))
plt.tight_layout()
fig.savefig(os.path.join(path_output_files, 'analysis_plots.png'), dpi=150)
plt.close(fig)

del data_mintflow, model, trainer, dfs, df_toinspect, df_filtered
gc.collect()
print(f"\nCompleted training for {run_label}")
