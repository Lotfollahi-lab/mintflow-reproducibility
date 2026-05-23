#!/usr/bin/env python

import gc
import gzip
import os
import pickle
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import mintflow
import mintflow.interface.base_interface
import mintflow.interface.perturbation.module_gen_micsizefactor
import mintflow.interface.perturbation.module_gen_stdata
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
import torch
from anndata._io.specs import registry as anndata_registry
from anndata._io.specs.methods import H5Array
from anndata._io.specs.registry import IOSpec
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test


# # Setup

def get_run_index():
    raw_value = os.environ.get('MINTFLOW_RUN_INDEX') or os.environ.get('LSB_JOBINDEX') or '1'
    try:
        run_index = int(raw_value)
    except ValueError as exc:
        raise ValueError(f'Invalid run index: {raw_value!r}') from exc

    if run_index < 1:
        raise ValueError(f'Run index must be >= 1, got {run_index}')

    return run_index


BASE_WORKFLOW_ROOT = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Final/Perturbations/Control_workflow'
)
BASE_RANDOM_SEED = 0
RUN_INDEX = get_run_index()
WORKFLOW_ROOT = BASE_WORKFLOW_ROOT / f'run_{RUN_INDEX:02d}'
RANDOM_SEED = BASE_RANDOM_SEED + RUN_INDEX - 1
RUN_DOWNSTREAM_ANALYSIS = False

MODEL_CHECKPOINT = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/ICB_TCell_Training/outputs/model_epoch_5.pt'
)
MACROPHAGE_REFERENCE_PATH = Path('/nfs/team361/cl35/Macs/adata_macrophages.h5ad')
PATH_CROPPED = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_TLS_cropped_combined'
)
PATH_ANNDATA = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_raw_per_section_ICB_TC_relabel'
)
NUM_TRAINING_EPOCHS = 8
SURVIVAL_EXPR_PATH = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Survival_analysis/HiSeqV2'
)
SURVIVAL_CLIN_PATH = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Survival_analysis/survival-KIRC_survival.txt'
)

RELABEL_MAP = {
    'CD4+ T lymphocyte': 'CD4+ T lymphocyte ICB',
    'CD8+ T lymphocyte': 'CD8+ T lymphocyte ICB',
    'DN T lymphocyte': 'DN T lymphocyte ICB',
}
TARGET_TCELL_LABELS = tuple(RELABEL_MAP)

FULL_REPLACEMENT_MODE = 'full_replacement'
DOSE_25_MODE = 'dose_response_25pct'
DOSE_50_MODE = 'dose_response_50pct'
DOSE_75_MODE = 'dose_response_75pct'
RANDOM_RELABEL_MODE = 'random_labels_weighted'

MODE_ORDER = [
    FULL_REPLACEMENT_MODE,
    DOSE_25_MODE,
    DOSE_50_MODE,
    DOSE_75_MODE,
    RANDOM_RELABEL_MODE,
]
MODE_CONFIGS = {
    FULL_REPLACEMENT_MODE: {
        'perturbation_mode': 'full_replacement',
        'dose_fraction': None,
        'random_label_sampling': 'weighted',
        'perturbation_label': 'Full replacement',
        'color': '#D55E00',
    },
    DOSE_25_MODE: {
        'perturbation_mode': 'dose_response',
        'dose_fraction': 0.25,
        'random_label_sampling': 'weighted',
        'perturbation_label': 'Dose 25%',
        'color': '#B79F00',
    },
    DOSE_50_MODE: {
        'perturbation_mode': 'dose_response',
        'dose_fraction': 0.5,
        'random_label_sampling': 'weighted',
        'perturbation_label': 'Dose 50%',
        'color': '#E69F00',
    },
    DOSE_75_MODE: {
        'perturbation_mode': 'dose_response',
        'dose_fraction': 0.75,
        'random_label_sampling': 'weighted',
        'perturbation_label': 'Dose 75%',
        'color': '#009E73',
    },
    RANDOM_RELABEL_MODE: {
        'perturbation_mode': 'random_labels',
        'dose_fraction': None,
        'random_label_sampling': 'weighted',
        'perturbation_label': 'Random relabel',
        'color': '#A65628',
    },
}

UNTREATED_SECTIONS = ['CV1-KID-0-FT-2', 'CV6-KID-0-FT-1', 'CV7-KID-0-FT-2-s3', 'CV9-KID-0-FT-2']
TREATED_SECTIONS = ['DI10-KID-0-FT-3', 'DI13-KID-0-FT-1']
GENERATED_SCORE_SECTIONS = UNTREATED_SECTIONS

BATCH_INDEX_MAP = {
    'CV6-KID-0-FT-1': 23,
    'CV9-KID-0-FT-2': 37,
    'DI10-KID-0-FT-3': 43,
    'DI13-KID-0-FT-1': 48,
    'CV1-KID-0-FO-1': 0,
    'CV1-KID-0-FT-2': 2,
    'CV5-KID-0-FO-1': 19,
    'CV7-KID-0-FT-2-s3': 30,
}

OBSKEY_CELLTYPE = 'level_3_cell_type'
NUM_GENERATED_REALISATIONS = 100
KWARGS_NEIGHBOURHOOD_GRAPH = {
    'spatial_key': 'spatial',
    'library_key': None,
    'set_diag': False,
    'delaunay': False,
    'n_neighs': 10,
}


def safe_name(value):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(value)).strip('_')


def save_current_figure(path, dpi=300):
    fig = plt.gcf()
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved plot -> {path}')


def save_plot_object(plot_obj, path, dpi=300):
    if hasattr(plot_obj, 'savefig'):
        plot_obj.savefig(path, dpi=dpi)
    else:
        plt.gcf().savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close('all')
    print(f'Saved plot -> {path}')


def save_adata(adata, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    print(f'Saved adata -> {path}')


def save_table(df, path, sep='\t', index=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=index)
    print(f'Saved table -> {path}')


def collect_memory():
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


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
    files = sorted([f for f in os.listdir(path_anndata) if f.endswith('.h5ad')])
    if not files:
        raise FileNotFoundError(f'No .h5ad files found in: {path_anndata}')

    num_sections = len(files)
    config_data_train, config_data_evaluation, config_model, config_training = (
        mintflow.get_default_configurations(
            num_tissue_sections_training=num_sections,
            num_tissue_sections_evaluation=num_sections,
        )
    )

    for idx, filename in enumerate(files, start=1):
        key = f'anndata{idx}'
        filepath = os.path.join(path_anndata, filename)

        config_data_train['list_tissue'][key]['file'] = filepath
        config_data_train['list_tissue'][key]['obskey_cell_type'] = 'level_3_cell_type'
        config_data_train['list_tissue'][key]['obskey_sliceid_to_checkUnique'] = 'donor'
        config_data_train['list_tissue'][key]['obskey_x'] = 'x_centroid'
        config_data_train['list_tissue'][key]['obskey_y'] = 'y_centroid'
        config_data_train['list_tissue'][key]['obskey_biological_batch_key'] = 'batch'
        config_data_train['list_tissue'][key]['config_dataloader_train']['width_window'] = 1368
        config_data_train['list_tissue'][key]['config_neighbourhood_graph'] = {
            'n_neighs': 10,
            'set_diag': 'False',
            'delaunay': 'False',
        }

        config_data_evaluation['list_tissue'][key]['file'] = filepath
        config_data_evaluation['list_tissue'][key]['obskey_cell_type'] = 'level_3_cell_type'
        config_data_evaluation['list_tissue'][key]['obskey_sliceid_to_checkUnique'] = 'donor'
        config_data_evaluation['list_tissue'][key]['obskey_x'] = 'x_centroid'
        config_data_evaluation['list_tissue'][key]['obskey_y'] = 'y_centroid'
        config_data_evaluation['list_tissue'][key]['obskey_biological_batch_key'] = 'batch'
        config_data_evaluation['list_tissue'][key]['config_dataloader_test']['width_window'] = 1368
        config_data_evaluation['list_tissue'][key]['config_neighbourhood_graph'] = {
            'n_neighs': 10,
            'set_diag': 'False',
            'delaunay': 'False',
        }

    config_model['coef_xbarint2notbatchID_loss'] = 1.0
    config_model['coef_xbarspl2notbatchID_loss'] = 1.0
    config_model['coef_flowmatchingloss'] = 0.0
    config_model['dict_qname_to_scaleandunweighted'] = (
        'impanddisentgl_int#0.1#True&impanddisentgl_spl#0.0#True'
        '&varphi_enc_int#0.0#True&varphi_enc_spl#0.0#True'
        '&z#0.0#True&sin#0.0#True&sout#0.0#True'
    )
    config_model['coef_loss_CTpredfromZ'] = 100

    config_training['num_training_epochs'] = num_training_epochs
    config_training['flag_use_GPU'] = 'True'
    config_training['flag_enable_wandb'] = 'False'
    config_training['annealing_decoder_XintXspl_coef_max'] = 0.01

    config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train)
    config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
    config_model = mintflow.verify_and_postprocess_config_model(
        config_model, num_tissue_sections=len(config_data_train)
    )
    config_training = mintflow.verify_and_postprocess_config_training(config_training)

    dict_all4_configs = {
        'config_data_train': config_data_train,
        'config_data_evaluation': config_data_evaluation,
        'config_model': config_model,
        'config_training': config_training,
    }
    _mintflow_bool_to_str(dict_all4_configs)
    data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
    return data_mintflow, dict_all4_configs


@anndata_registry._REGISTRY.register_read(H5Array, IOSpec('null', '0.1.0'))
def _read_h5ad_null_compat(elem, *, _reader=None):
    return None


def unique_preserve_order(values):
    seen = set()
    ordered = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def assign_level3_labels(adata, new_labels):
    new_series = pd.Series(new_labels, index=adata.obs_names, dtype='object').astype(str)
    old_series = adata.obs['level_3_cell_type']

    if pd.api.types.is_categorical_dtype(old_series):
        categories = unique_preserve_order(old_series.astype(str).tolist() + new_series.tolist())
        adata.obs['level_3_cell_type'] = pd.Series(
            pd.Categorical(new_series, categories=categories),
            index=adata.obs_names,
        ).cat.remove_unused_categories()
    else:
        adata.obs['level_3_cell_type'] = new_series.to_numpy()

    return adata


def is_tcell_like_label(label):
    label = str(label)
    lower_label = label.lower()
    return (
        label in RELABEL_MAP
        or label in RELABEL_MAP.values()
        or 't lymphocyte' in lower_label
        or 't cell' in lower_label
        or 't-cell' in lower_label
    )


def generate_perturbed_level3_labels(section_name, adata, rng, mode_config):
    old_labels = pd.Series(adata.obs['level_3_cell_type'], index=adata.obs_names).astype(str)
    target_mask = old_labels.isin(TARGET_TCELL_LABELS)
    new_labels = old_labels.copy()

    summary = {
        'section': section_name,
        'perturbation_mode': mode_config['perturbation_mode'],
        'dose_fraction': (
            mode_config['dose_fraction']
            if mode_config['perturbation_mode'] == 'dose_response'
            else np.nan
        ),
        'random_label_sampling': (
            mode_config['random_label_sampling']
            if mode_config['perturbation_mode'] == 'random_labels'
            else ''
        ),
        'n_section_cells': int(old_labels.shape[0]),
        'n_target_t_cells': int(target_mask.sum()),
    }

    if mode_config['perturbation_mode'] == 'full_replacement':
        new_labels.loc[target_mask] = new_labels.loc[target_mask].map(RELABEL_MAP)

    elif mode_config['perturbation_mode'] == 'dose_response':
        for source_label, target_label in RELABEL_MAP.items():
            label_indices = np.flatnonzero(old_labels.to_numpy() == source_label)
            ordered_indices = rng.permutation(label_indices)
            n_selected = min(
                len(ordered_indices),
                int(np.rint(mode_config['dose_fraction'] * len(ordered_indices))),
            )
            selected_indices = ordered_indices[:n_selected]
            new_labels.iloc[selected_indices] = target_label

    elif mode_config['perturbation_mode'] == 'random_labels':
        eligible_counts = old_labels.loc[~old_labels.map(is_tcell_like_label)].value_counts()
        if eligible_counts.empty and target_mask.any():
            raise ValueError(
                f'No eligible non-T level_3 labels were found in section {section_name} '
                'for random-label assignment.'
            )

        if target_mask.any():
            eligible_labels = eligible_counts.index.to_numpy(dtype=object)
            sampled_labels = rng.choice(
                eligible_labels,
                size=int(target_mask.sum()),
                replace=True,
                p=(
                    eligible_counts.to_numpy(dtype=float) / eligible_counts.sum()
                    if mode_config['random_label_sampling'] == 'weighted'
                    else None
                ),
            )
            new_labels.loc[target_mask] = sampled_labels

        summary['n_random_label_candidates'] = int(eligible_counts.shape[0])

    else:
        raise ValueError(f"Unsupported perturbation mode: {mode_config['perturbation_mode']}")

    changed_mask = new_labels != old_labels
    summary['n_changed_cells'] = int(changed_mask.sum())

    change_counts = (
        pd.DataFrame({'old_label': old_labels, 'new_label': new_labels})
        .loc[changed_mask]
        .value_counts(sort=True)
        .rename('n_cells')
        .reset_index()
    )
    if change_counts.empty:
        change_counts = pd.DataFrame(columns=['old_label', 'new_label', 'n_cells'])

    change_counts.insert(0, 'section', section_name)
    change_counts.insert(1, 'perturbation_mode', mode_config['perturbation_mode'])

    return new_labels, summary, change_counts


def mean_generated_expression(result_gen):
    expr_sum = None
    n_realisations = 0
    for realisation in result_gen['list_generated_realisations_ie_expressions']:
        generated = np.asarray(realisation['MintFLow_Generated_Xmic'])
        if expr_sum is None:
            expr_sum = generated.astype(np.float64, copy=True)
        else:
            expr_sum += generated
        n_realisations += 1

    if n_realisations == 0:
        raise ValueError('No generated realisations were returned by mintflow.generate_insilico_ST_data')

    return expr_sum / n_realisations


def build_mcc_mask_table(adata, section_name, condition_name, changed_mask):
    if 'cell_id' in adata.obs.columns:
        cell_ids = adata.obs['cell_id'].astype(str).to_numpy()
    else:
        cell_ids = adata.obs_names.astype(str).to_numpy()

    df = pd.DataFrame(
        {
            'section': section_name,
            'condition': condition_name,
            'obs_name': adata.obs_names.astype(str),
            'cell_id': cell_ids,
            'mcc_changed': np.asarray(changed_mask, dtype=bool),
        }
    )

    for col in ['batch', 'level_3_cell_type', 'level_4_cell_type']:
        if col in adata.obs.columns:
            df[col] = adata.obs[col].to_numpy()

    if 'spatial' in adata.obsm:
        spatial = np.asarray(adata.obsm['spatial'])
        if spatial.ndim == 2 and spatial.shape[1] >= 1:
            df['spatial_x'] = spatial[:, 0]
        if spatial.ndim == 2 and spatial.shape[1] >= 2:
            df['spatial_y'] = spatial[:, 1]

    return df


def build_generated_expression_adata(adata):
    if 'MintFLow_Generated_Xmic' not in adata.obsm:
        raise KeyError("Missing `obsm['MintFLow_Generated_Xmic']` in saved adata.")

    score_adata = sc.AnnData(
        X=np.asarray(adata.obsm['MintFLow_Generated_Xmic'], dtype=np.float64),
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    score_adata.obs_names = adata.obs_names.copy()
    score_adata.var_names = adata.var_names.copy()
    score_adata.var_names_make_unique()
    return score_adata


# # Scanpy and device

sc.settings.verbosity = 3
sc.settings.seed = RANDOM_SEED
print('GPU Available:', torch.cuda.is_available())
print(f'Run index: {RUN_INDEX}')
print(f'Random seed: {RANDOM_SEED}')

device = torch.device('cpu')
if torch.cuda.is_available():
    print('GPU Name:', torch.cuda.get_device_name(0))
    device = torch.device('cuda')


# # Output layout

workflow_root = WORKFLOW_ROOT
workflow_root.mkdir(parents=True, exist_ok=True)
print(f'Workflow output root: {workflow_root}')

shared_root = workflow_root / 'shared_unperturbed'
shared_label_dir = shared_root / 'label_adatas'
shared_generated_dir = shared_root / 'generated_adatas'
observed_deg_output_dir = workflow_root / 'shared_observed_degs'

mode_dirs = {}
for mode_name in MODE_ORDER:
    mode_root = workflow_root / mode_name
    mode_dirs[mode_name] = {
        'root': mode_root,
        'label_adatas': mode_root / 'label_adatas',
        'generated_adatas': mode_root / 'generated_adatas',
        'mcc_masks': mode_root / 'mcc_masks',
        'deg_tables': mode_root / 'deg_tables',
        'plots': mode_root / 'plots',
        'signature_scores': mode_root / 'observed_deg_signature_scores',
        'results_pickle': mode_root / 'results.pkl.gz',
    }

shared_label_dir.mkdir(parents=True, exist_ok=True)
shared_generated_dir.mkdir(parents=True, exist_ok=True)
observed_deg_output_dir.mkdir(parents=True, exist_ok=True)


# # Load checkpoint

model = torch.load(MODEL_CHECKPOINT, map_location='cpu', weights_only=False)
model.to(device)
print(f'Loaded checkpoint: {MODEL_CHECKPOINT}')

data_mintflow, dict_all4_configs = build_mintflow_context(
    PATH_ANNDATA, num_training_epochs=NUM_TRAINING_EPOCHS
)
print(f'Built MintFlow context from: {PATH_ANNDATA}')


# # Load reference and sections

adata_macrophages = sc.read_h5ad(MACROPHAGE_REFERENCE_PATH)
ref_df = adata_macrophages.obs[['batch', 'Annotation_TAMs_merged']].copy()
ref_df['cell_id'] = adata_macrophages.obs_names.astype(str)
ref_level4 = ref_df.set_index(['batch', 'cell_id'])['Annotation_TAMs_merged']

n_dup = ref_level4.index.duplicated(keep='first').sum()
if n_dup > 0:
    print(f'Found {n_dup} duplicated (batch, cell_id) entries in reference; keeping first.')
ref_level4 = ref_level4[~ref_level4.index.duplicated(keep='first')]

del ref_df
del adata_macrophages
collect_memory()

files = sorted([f for f in os.listdir(PATH_CROPPED) if f.endswith('.h5ad')])
section_names = [f.replace('_crop.h5ad', '') for f in files]
print('Samples:', section_names)

base_adatas = {}
for section_name in section_names:
    adata_sample = sc.read_h5ad(PATH_CROPPED / f'{section_name}_crop.h5ad')

    if 'batch' in adata_sample.obs.columns:
        batch_vals = adata_sample.obs['batch'].astype(str).to_numpy()
    else:
        batch_vals = np.repeat(section_name, adata_sample.n_obs)

    key = pd.MultiIndex.from_arrays(
        [batch_vals, adata_sample.obs_names.astype(str)],
        names=['batch', 'cell_id'],
    )
    adata_sample.obs['level_4_cell_type'] = pd.Categorical(ref_level4.reindex(key).to_numpy())
    adata_sample.obs['section'] = section_name
    base_adatas[section_name] = adata_sample

    n_missing = int(adata_sample.obs['level_4_cell_type'].isna().sum())
    print(f'Loaded {section_name}: {adata_sample.n_obs} cells (unmatched labels: {n_missing})')


# # Prepare label-stage adatas

shared_rows = []
for section_name, adata in base_adatas.items():
    adata_original = adata.copy()
    save_adata(adata_original, shared_label_dir / f'{safe_name(section_name)}__unperturbed.h5ad')
    shared_rows.append({'section': section_name, 'n_section_cells': int(adata_original.n_obs)})

save_table(pd.DataFrame(shared_rows), shared_root / 'section_summary.tsv')


def prepare_mode_label_adatas(mode_name):
    mode_dir = mode_dirs[mode_name]
    mode_dir['label_adatas'].mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(RANDOM_SEED)
    perturbation_section_summaries = []
    perturbation_label_change_tables = []

    for section_name, adata in base_adatas.items():
        adata_perturbed = adata.copy()

        new_labels, section_summary, change_counts = generate_perturbed_level3_labels(
            section_name=section_name,
            adata=adata_perturbed,
            rng=rng,
            mode_config=MODE_CONFIGS[mode_name],
        )
        assign_level3_labels(adata_perturbed, new_labels)

        save_adata(
            adata_perturbed,
            mode_dir['label_adatas'] / f'{safe_name(section_name)}__perturbed.h5ad',
        )

        perturbation_section_summaries.append(section_summary)
        perturbation_label_change_tables.append(change_counts)

    perturbation_section_summary_df = pd.DataFrame(perturbation_section_summaries)
    if perturbation_label_change_tables:
        perturbation_label_change_df = pd.concat(perturbation_label_change_tables, ignore_index=True)
    else:
        perturbation_label_change_df = pd.DataFrame(
            columns=['section', 'perturbation_mode', 'old_label', 'new_label', 'n_cells']
        )

    save_table(perturbation_section_summary_df, mode_dir['root'] / 'perturbation_section_summary.tsv')
    save_table(perturbation_label_change_df, mode_dir['root'] / 'perturbation_label_changes.tsv.gz')
    return mode_name


with ThreadPoolExecutor(max_workers=6) as executor:
    futures = [executor.submit(prepare_mode_label_adatas, mode_name) for mode_name in MODE_ORDER]
    for future in futures:
        print(f'Prepared label-stage adatas for mode: {future.result()}')

del base_adatas
collect_memory()


# # Unperturbed generation

for section_name in section_names:
    print(f"\n{'=' * 60}")
    print(f'Processing shared unperturbed baseline: {section_name}')
    print(f"{'=' * 60}")

    adata_orig = sc.read_h5ad(shared_label_dir / f'{safe_name(section_name)}__unperturbed.h5ad')
    batch_idx = BATCH_INDEX_MAP[section_name]
    adata_output_path = shared_generated_dir / f'{safe_name(section_name)}__unperturbed.h5ad'

    adata_orig.uns = {}
    sq.gr.spatial_neighbors(adata=adata_orig, **KWARGS_NEIGHBOURHOOD_GRAPH)

    print(f'  Generating expression for shared unperturbed tissue (batch_index={batch_idx})...')
    result_gen_orig = mintflow.generate_insilico_ST_data(
        adata=adata_orig,
        obskey_celltype=OBSKEY_CELLTYPE,
        obspkey_neighbourhood_graph='spatial_connectivities',
        device=device,
        batch_index_trainingdata=batch_idx,
        num_generated_realisations=NUM_GENERATED_REALISATIONS,
        model=model,
        data_mintflow=data_mintflow,
        dict_all4_configs=dict_all4_configs,
        estimate_spatial_sizefactors_on_sections=[batch_idx],
    )

    adata_orig.obsm['MintFLow_Generated_Xmic'] = mean_generated_expression(result_gen_orig)
    adata_orig.obsm['MintFLow_np_MCC'] = np.asarray(result_gen_orig['np_MCC'])
    save_adata(adata_orig, adata_output_path)

    del adata_orig
    del result_gen_orig
    collect_memory()


# # Perturbed and control generation

for mode_name in MODE_ORDER:
    print(f"\n{'#' * 80}")
    print(f'Generating mode: {mode_name}')
    print(f"{'#' * 80}")

    mode_dir = mode_dirs[mode_name]
    mode_dir['generated_adatas'].mkdir(parents=True, exist_ok=True)
    mode_dir['mcc_masks'].mkdir(parents=True, exist_ok=True)
    results = {}

    for section_name in section_names:
        print(f"\n{'=' * 60}")
        print(f'Processing: {section_name} ({mode_name})')
        print(f"{'=' * 60}")

        adata_orig = sc.read_h5ad(shared_generated_dir / f'{safe_name(section_name)}__unperturbed.h5ad')
        adata_pert = sc.read_h5ad(mode_dir['label_adatas'] / f'{safe_name(section_name)}__perturbed.h5ad')
        batch_idx = BATCH_INDEX_MAP[section_name]

        if 'MintFLow_np_MCC' not in adata_orig.obsm:
            raise KeyError(f'Missing shared baseline MCC array for {section_name}')

        list_cell_id_orig = adata_orig.obs['cell_id'].astype(str).tolist()
        np_mcc_orig = np.asarray(adata_orig.obsm['MintFLow_np_MCC'])
        adata_original_path = shared_generated_dir / f'{safe_name(section_name)}__unperturbed.h5ad'
        adata_perturbed_path = mode_dir['generated_adatas'] / f'{safe_name(section_name)}__perturbed.h5ad'

        adata_pert.uns = {}
        sq.gr.spatial_neighbors(adata=adata_pert, **KWARGS_NEIGHBOURHOOD_GRAPH)

        print('  Generating expression for perturbed tissue...')
        result_gen_pert = mintflow.generate_insilico_ST_data(
            adata=adata_pert,
            obskey_celltype=OBSKEY_CELLTYPE,
            obspkey_neighbourhood_graph='spatial_connectivities',
            device=device,
            batch_index_trainingdata=batch_idx,
            num_generated_realisations=NUM_GENERATED_REALISATIONS,
            model=model,
            data_mintflow=data_mintflow,
            dict_all4_configs=dict_all4_configs,
            estimate_spatial_sizefactors_on_sections=[batch_idx],
        )

        adata_pert.obsm['MintFLow_Generated_Xmic'] = mean_generated_expression(result_gen_pert)
        adata_pert.obsm['MintFLow_np_MCC'] = np.asarray(result_gen_pert['np_MCC'])
        list_cell_id_pert = adata_pert.obs['cell_id'].astype(str).tolist()
        np_mcc_pert = result_gen_pert['np_MCC'].copy()
        save_adata(adata_pert, adata_perturbed_path)

        dict_cell_id_to_idx_orig = {cid: idx for idx, cid in enumerate(list_cell_id_orig)}
        dict_cell_id_to_idx_pert = {cid: idx for idx, cid in enumerate(list_cell_id_pert)}

        set_orig = set(list_cell_id_orig)
        set_pert = set(list_cell_id_pert)
        all_cells = set_orig | set_pert
        common_cells = set_orig & set_pert

        dict_cell_id_to_mcc_changed = {}
        for cid in all_cells:
            if cid not in common_cells:
                dict_cell_id_to_mcc_changed[cid] = False
            else:
                dict_cell_id_to_mcc_changed[cid] = not np.allclose(
                    np_mcc_orig[dict_cell_id_to_idx_orig[cid]],
                    np_mcc_pert[dict_cell_id_to_idx_pert[cid]],
                )

        print(
            'For {} percent of cells, the MCC has changed due to perturbation.'.format(
                np.round(100.0 * np.mean(list(dict_cell_id_to_mcc_changed.values())), 3)
            )
        )

        changed_orig = np.array(
            [dict_cell_id_to_mcc_changed.get(cid, False) for cid in list_cell_id_orig],
            dtype=bool,
        )
        changed_pert = np.array(
            [dict_cell_id_to_mcc_changed.get(cid, False) for cid in list_cell_id_pert],
            dtype=bool,
        )

        mcc_mask_df = pd.concat(
            [
                build_mcc_mask_table(adata_orig, section_name, 'unperturbed', changed_orig),
                build_mcc_mask_table(adata_pert, section_name, 'perturbed', changed_pert),
            ],
            ignore_index=True,
        )
        mcc_mask_path = mode_dir['mcc_masks'] / f'{safe_name(section_name)}__mcc_masks.tsv.gz'
        save_table(mcc_mask_df, mcc_mask_path)

        results[section_name] = {
            'adata_original_path': str(adata_original_path),
            'adata_perturbed_path': str(adata_perturbed_path),
            'mcc_mask_path': str(mcc_mask_path),
            'dict_cellID_to_MCCchanged': dict_cell_id_to_mcc_changed,
        }

        del adata_orig
        del adata_pert
        del result_gen_pert
        del np_mcc_pert
        collect_memory()

    with gzip.open(mode_dir['results_pickle'], 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved results pickle -> {mode_dir['results_pickle']}")
