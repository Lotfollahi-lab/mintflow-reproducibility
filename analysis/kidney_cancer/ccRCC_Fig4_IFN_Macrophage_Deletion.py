#!/usr/bin/env python

import gc
import gzip
import math
import os
import pickle
import re
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
from matplotlib.lines import Line2D


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
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Final/Perturbations/Control_workflow_mac_deletion'
)
BASE_RANDOM_SEED = 0
RUN_INDEX = get_run_index()
WORKFLOW_ROOT = BASE_WORKFLOW_ROOT / f'run_{RUN_INDEX:02d}'
RANDOM_SEED = BASE_RANDOM_SEED + RUN_INDEX - 1

MODEL_CHECKPOINT = Path(
    '/nfs/team361/am74/mintflow/revision/train_model_all_data/window_sz_2280_alldatawithbatcheffect/8epochs_save_h5ad/new_hyperparametes/run2/outputs/model_epoch_7.pt'
)
MACROPHAGE_REFERENCE_PATH = Path('/nfs/team361/cl35/Macs/adata_macrophages.h5ad')
TCELL_REFERENCE_PATH = Path('/nfs/team361/cl35/Tcells/adata_lymphocyte_annotated.h5ad')
PATH_CROPPED = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_TLS_cropped_combined'
)
PATH_ANNDATA = Path(
    '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_raw_per_section'
)
NUM_TRAINING_EPOCHS = 8

TARGET_SECTIONS = [
    'CV1-KID-0-FT-2',
    'CV6-KID-0-FT-1',
    'CV7-KID-0-FT-2-s3',
    'CV9-KID-0-FT-2',
]

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

IFN_TARGET_LABEL = 'IFN+'
OBSKEY_CELLTYPE = 'level_3_cell_type'
NUM_GENERATED_REALISATIONS = 100
KWARGS_NEIGHBOURHOOD_GRAPH = {
    'spatial_key': 'spatial',
    'library_key': None,
    'set_diag': False,
    'delaunay': False,
    'n_neighs': 10,
}


UNPERTURBED_MODE = 'unperturbed'
FULL_DELETION_MODE = 'full_deletion'
DOSE_25_MODE = 'dose_response_25pct'
DOSE_50_MODE = 'dose_response_50pct'
DOSE_75_MODE = 'dose_response_75pct'
RANDOM_DELETION_MODE = 'random_deletion'

MODE_ORDER = [
    FULL_DELETION_MODE,
    DOSE_25_MODE,
    DOSE_50_MODE,
    DOSE_75_MODE,
    RANDOM_DELETION_MODE,
]

MODE_CONFIGS = {
    FULL_DELETION_MODE: {
        'perturbation_mode': 'full_deletion',
        'dose_fraction': None,
        'perturbation_label': 'Full deletion',
        'color': '#D55E00',
    },
    DOSE_25_MODE: {
        'perturbation_mode': 'dose_response',
        'dose_fraction': 0.25,
        'perturbation_label': 'Dose 25%',
        'color': '#B79F00',
    },
    DOSE_50_MODE: {
        'perturbation_mode': 'dose_response',
        'dose_fraction': 0.50,
        'perturbation_label': 'Dose 50%',
        'color': '#E69F00',
    },
    DOSE_75_MODE: {
        'perturbation_mode': 'dose_response',
        'dose_fraction': 0.75,
        'perturbation_label': 'Dose 75%',
        'color': '#009E73',
    },
    RANDOM_DELETION_MODE: {
        'perturbation_mode': 'random_deletion',
        'dose_fraction': None,
        'perturbation_label': 'Random deletion',
        'color': '#A65628',
    },
}

PLOT_CONDITION_ORDER = [
    'Unperturbed',
    'Dose 25%',
    'Dose 50%',
    'Dose 75%',
    'Full deletion',
    'Random deletion',
]


def safe_name(value):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(value)).strip('_')


def save_current_figure(path, dpi=300):
    fig = plt.gcf()
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
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


def ensure_cell_id_and_batch(adata, section_name=None):
    if 'cell_id' not in adata.obs.columns:
        if 'Cell_ID' in adata.obs.columns:
            adata.obs['cell_id'] = adata.obs['Cell_ID'].astype(str)
        else:
            adata.obs['cell_id'] = adata.obs_names.astype(str)

    if 'batch' not in adata.obs.columns:
        if section_name is None:
            raise KeyError('`batch` not found in adata.obs')
        adata.obs['batch'] = section_name

    adata.obs['cell_id'] = adata.obs['cell_id'].astype(str).str.strip()
    adata.obs['batch'] = adata.obs['batch'].astype(str).str.strip()
    return adata


def mean_generated_expression(result_gen, expression_key='MintFLow_Generated_Xmic'):
    expr_sum = None
    n_realisations = 0
    for realisation in result_gen['list_generated_realisations_ie_expressions']:
        if expression_key not in realisation:
            raise KeyError(f'Missing generated expression key {expression_key!r} in MintFlow output')
        generated = np.asarray(realisation[expression_key])
        if expr_sum is None:
            expr_sum = generated.astype(np.float64, copy=True)
        else:
            expr_sum += generated
        n_realisations += 1

    if n_realisations == 0:
        raise ValueError('No generated realisations were returned by mintflow.generate_insilico_ST_data')

    return expr_sum / n_realisations


def build_generated_expression_adata(adata):
    for key in ('MintFlow_Generated_Xint', 'MintFLow_Generated_Xmic'):
        if key not in adata.obsm:
            raise KeyError(f"Missing `obsm[{key!r}]` in saved adata.")

    # total generated expression = intrinsic (Xint) + microenvironment (Xmic)
    score_adata = sc.AnnData(
        X=(
            np.asarray(adata.obsm['MintFlow_Generated_Xint'], dtype=np.float64)
            + np.asarray(adata.obsm['MintFLow_Generated_Xmic'], dtype=np.float64)
        ),
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    score_adata.obs_names = adata.obs_names.copy()
    score_adata.var_names = adata.var_names.copy()
    score_adata.var_names_make_unique()
    return score_adata


def match_signature_genes(input_genes, var_names):
    lookup = {}
    for gene in var_names:
        lookup.setdefault(str(gene).upper(), str(gene))

    matched = []
    missing = []
    for gene in input_genes:
        match = lookup.get(str(gene).upper())
        if match is None:
            missing.append(gene)
        else:
            matched.append(match)

    return unique_preserve_order(matched), missing


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


def load_tcell_level4_reference():
    adata_tcells = sc.read_h5ad(TCELL_REFERENCE_PATH)
    adata_tcells = ensure_cell_id_and_batch(adata_tcells)
    ref_df = adata_tcells.obs[['batch', 'cell_id', 'Annotation_XMic_Tcells']].copy()
    ref_df['batch'] = ref_df['batch'].astype(str).str.strip()
    ref_df['cell_id'] = ref_df['cell_id'].astype(str).str.strip()
    ref_level4 = ref_df.set_index(['batch', 'cell_id'])['Annotation_XMic_Tcells']

    n_dup = ref_level4.index.duplicated(keep='first').sum()
    if n_dup > 0:
        print(f'Found {n_dup} duplicated (batch, cell_id) entries in T-cell reference; keeping first.')
    ref_level4 = ref_level4[~ref_level4.index.duplicated(keep='first')]

    del ref_df
    del adata_tcells
    collect_memory()
    return ref_level4


def load_ifn_macrophage_targets():
    adata_macs = sc.read_h5ad(MACROPHAGE_REFERENCE_PATH)
    adata_macs = ensure_cell_id_and_batch(adata_macs)

    target_df = adata_macs.obs.loc[
        adata_macs.obs['Annotation_TAMs_merged'].astype(str) == IFN_TARGET_LABEL,
        ['batch', 'cell_id'],
    ].copy()
    target_df['batch'] = target_df['batch'].astype(str).str.strip()
    target_df['cell_id'] = target_df['cell_id'].astype(str).str.strip()
    target_df = target_df.drop_duplicates().sort_values(['batch', 'cell_id']).reset_index(drop=True)
    target_keys = pd.MultiIndex.from_frame(target_df[['batch', 'cell_id']])

    del adata_macs
    collect_memory()
    return target_keys, target_df


def load_base_adatas(ref_level4):
    base_adatas = {}

    for section_name in TARGET_SECTIONS:
        path = PATH_CROPPED / f'{section_name}_crop.h5ad'
        if not path.exists():
            raise FileNotFoundError(f'Missing section file: {path}')

        adata = sc.read_h5ad(path)
        adata = ensure_cell_id_and_batch(adata, section_name=section_name)

        key = pd.MultiIndex.from_arrays(
            [adata.obs['batch'].astype(str).to_numpy(), adata.obs['cell_id'].astype(str).to_numpy()],
            names=['batch', 'cell_id'],
        )
        adata.obs['level_4_cell_type'] = pd.Categorical(ref_level4.reindex(key).to_numpy())
        adata.obs['section'] = section_name

        base_adatas[section_name] = adata
        n_missing = int(adata.obs['level_4_cell_type'].isna().sum())
        print(f'Loaded {section_name}: {adata.n_obs} cells (unmatched level_4 labels: {n_missing})')

    return base_adatas


def build_ifn_target_mask(adata, target_keys):
    obs_keys = pd.MultiIndex.from_frame(adata.obs[['batch', 'cell_id']])
    return pd.Series(obs_keys.isin(target_keys), index=adata.obs_names)


def build_section_perturbation_plans(base_adatas, target_keys):
    rng = np.random.default_rng(RANDOM_SEED)
    plan_summary_rows = []
    target_rows = []
    plans = {}

    for section_name in TARGET_SECTIONS:
        adata = base_adatas[section_name]
        target_mask = build_ifn_target_mask(adata, target_keys)

        target_obs_names = adata.obs_names[target_mask.to_numpy()].astype(str).to_numpy()
        ordered_target_obs_names = rng.permutation(target_obs_names)
        ordered_random_delete_obs_names = rng.permutation(adata.obs_names.astype(str).to_numpy())

        retained_mask = ~target_mask.to_numpy()
        retained_obs_names = adata.obs_names[retained_mask].astype(str).to_numpy()

        target_df = adata.obs.loc[adata.obs_names[target_mask.to_numpy()], ['batch', 'cell_id']].copy()
        target_df.insert(0, 'section', section_name)
        target_df.insert(1, 'obs_name', target_df.index.astype(str))
        target_rows.append(target_df.reset_index(drop=True))

        plans[section_name] = {
            'ordered_target_obs_names': ordered_target_obs_names,
            'ordered_random_delete_obs_names': ordered_random_delete_obs_names,
            'retained_obs_names_after_full_deletion': retained_obs_names,
            'target_obs_name_set': set(target_obs_names.tolist()),
            'n_target_cells': int(target_mask.sum()),
            'n_section_cells': int(adata.n_obs),
        }

        plan_summary_rows.append(
            {
                'section': section_name,
                'n_section_cells': int(adata.n_obs),
                'n_target_ifn_macrophages': int(target_mask.sum()),
            }
        )

    target_cells_df = (
        pd.concat(target_rows, ignore_index=True)
        if target_rows
        else pd.DataFrame(columns=['section', 'obs_name', 'batch', 'cell_id'])
    )
    return plans, pd.DataFrame(plan_summary_rows), target_cells_df


def determine_deleted_obs_names(mode_name, plan):
    if mode_name == FULL_DELETION_MODE:
        return plan['ordered_target_obs_names']

    if mode_name in {DOSE_25_MODE, DOSE_50_MODE, DOSE_75_MODE}:
        dose_fraction = MODE_CONFIGS[mode_name]['dose_fraction']
        n_selected = min(
            plan['n_target_cells'],
            int(np.rint(dose_fraction * plan['n_target_cells'])),
        )
        return plan['ordered_target_obs_names'][:n_selected]

    if mode_name == RANDOM_DELETION_MODE:
        return plan['ordered_random_delete_obs_names'][:plan['n_target_cells']]

    raise ValueError(f'Unsupported mode: {mode_name}')


def build_deleted_cells_table(adata_before, section_name, mode_name, deleted_obs_names, plan):
    deleted_obs_name_set = set(np.asarray(deleted_obs_names, dtype=object).tolist())
    deleted_mask = adata_before.obs_names.astype(str).isin(deleted_obs_name_set)
    if not np.any(deleted_mask):
        return pd.DataFrame(
            columns=[
                'section',
                'mode_name',
                'perturbation_label',
                'obs_name',
                'cell_id',
                'batch',
                'level_3_cell_type',
                'level_4_cell_type',
                'is_target_ifn_macrophage',
            ]
        )

    deleted_df = adata_before.obs.loc[adata_before.obs_names[deleted_mask]].copy()
    deleted_df['section'] = section_name
    deleted_df['mode_name'] = mode_name
    deleted_df['perturbation_label'] = MODE_CONFIGS[mode_name]['perturbation_label']
    deleted_df['obs_name'] = deleted_df.index.astype(str)
    deleted_df['is_target_ifn_macrophage'] = deleted_df['obs_name'].isin(plan['target_obs_name_set'])

    keep_cols = [
        'section',
        'mode_name',
        'perturbation_label',
        'obs_name',
        'cell_id',
        'batch',
        'level_3_cell_type',
        'level_4_cell_type',
        'is_target_ifn_macrophage',
    ]
    return deleted_df.loc[:, keep_cols].reset_index(drop=True)


def prepare_mode_input_adatas(mode_name, base_adatas, section_plans, mode_dir):
    input_dir = mode_dir['input_adatas']
    input_dir.mkdir(parents=True, exist_ok=True)

    section_rows = []
    deleted_tables = []

    for section_name in TARGET_SECTIONS:
        adata_before = base_adatas[section_name]
        plan = section_plans[section_name]
        deleted_obs_names = determine_deleted_obs_names(mode_name, plan)
        deleted_obs_name_set = set(np.asarray(deleted_obs_names, dtype=object).tolist())
        delete_mask = adata_before.obs_names.astype(str).isin(deleted_obs_name_set)

        adata_after = adata_before[~delete_mask].copy()

        save_adata(
            adata_after,
            input_dir / f'{safe_name(section_name)}__perturbed.h5ad',
        )

        deleted_table = build_deleted_cells_table(
            adata_before=adata_before,
            section_name=section_name,
            mode_name=mode_name,
            deleted_obs_names=deleted_obs_names,
            plan=plan,
        )
        deleted_tables.append(deleted_table)

        n_deleted_target = int(deleted_table['is_target_ifn_macrophage'].sum()) if not deleted_table.empty else 0
        section_rows.append(
            {
                'section': section_name,
                'mode_name': mode_name,
                'perturbation_mode': MODE_CONFIGS[mode_name]['perturbation_mode'],
                'perturbation_label': MODE_CONFIGS[mode_name]['perturbation_label'],
                'dose_fraction': (
                    MODE_CONFIGS[mode_name]['dose_fraction']
                    if MODE_CONFIGS[mode_name]['perturbation_mode'] == 'dose_response'
                    else np.nan
                ),
                'n_section_cells': plan['n_section_cells'],
                'n_target_ifn_macrophages': plan['n_target_cells'],
                'n_deleted_cells': int(delete_mask.sum()),
                'n_deleted_target_ifn_macrophages': n_deleted_target,
                'n_remaining_cells': int(adata_after.n_obs),
            }
        )

        del adata_after
        collect_memory()

    section_summary_df = pd.DataFrame(section_rows)
    deleted_cells_df = (
        pd.concat(deleted_tables, ignore_index=True)
        if deleted_tables
        else pd.DataFrame()
    )

    save_table(section_summary_df, mode_dir['root'] / 'perturbation_section_summary.tsv')
    save_table(deleted_cells_df, mode_dir['root'] / 'deleted_cells.tsv.gz')


def generate_shared_unperturbed(section_names, shared_input_dir, shared_generated_dir, model, data_mintflow, dict_all4_configs, device):
    for section_name in section_names:
        print(f"\n{'=' * 60}")
        print(f'Processing shared unperturbed baseline: {section_name}')
        print(f"{'=' * 60}")

        adata_orig = sc.read_h5ad(shared_input_dir / f'{safe_name(section_name)}__unperturbed.h5ad')
        batch_idx = BATCH_INDEX_MAP[section_name]

        adata_orig.uns = {}
        sq.gr.spatial_neighbors(adata=adata_orig, **KWARGS_NEIGHBOURHOOD_GRAPH)

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

        adata_orig.obsm['MintFLow_Generated_Xmic'] = mean_generated_expression(
            result_gen_orig,
            expression_key='MintFLow_Generated_Xmic',
        )
        adata_orig.obsm['MintFlow_Generated_Xint'] = mean_generated_expression(
            result_gen_orig,
            expression_key='MintFlow_Generated_Xint',
        )
        adata_orig.obsm['MintFlow_Generated_Xint_plus_Xmic'] = (
            np.asarray(adata_orig.obsm['MintFlow_Generated_Xint'], dtype=np.float64)
            + np.asarray(adata_orig.obsm['MintFLow_Generated_Xmic'], dtype=np.float64)
        )
        adata_orig.obsm['MintFLow_np_MCC'] = np.asarray(result_gen_orig['np_MCC'])
        save_adata(
            adata_orig,
            shared_generated_dir / f'{safe_name(section_name)}__unperturbed.h5ad',
        )

        del adata_orig
        del result_gen_orig
        collect_memory()


def generate_mode_outputs(section_names, mode_name, mode_dir, shared_generated_dir, model, data_mintflow, dict_all4_configs, device):
    mode_dir['generated_adatas'].mkdir(parents=True, exist_ok=True)
    mode_dir['mcc_masks'].mkdir(parents=True, exist_ok=True)
    results = {}

    for section_name in section_names:
        print(f"\n{'=' * 60}")
        print(f'Processing: {section_name} ({mode_name})')
        print(f"{'=' * 60}")

        adata_orig = sc.read_h5ad(shared_generated_dir / f'{safe_name(section_name)}__unperturbed.h5ad')
        adata_pert = sc.read_h5ad(mode_dir['input_adatas'] / f'{safe_name(section_name)}__perturbed.h5ad')
        batch_idx = BATCH_INDEX_MAP[section_name]

        if 'MintFLow_np_MCC' not in adata_orig.obsm:
            raise KeyError(f'Missing shared baseline MCC array for {section_name}')

        adata_pert.uns = {}
        sq.gr.spatial_neighbors(adata=adata_pert, **KWARGS_NEIGHBOURHOOD_GRAPH)

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

        adata_pert.obsm['MintFLow_Generated_Xmic'] = mean_generated_expression(
            result_gen_pert,
            expression_key='MintFLow_Generated_Xmic',
        )
        adata_pert.obsm['MintFlow_Generated_Xint'] = mean_generated_expression(
            result_gen_pert,
            expression_key='MintFlow_Generated_Xint',
        )
        adata_pert.obsm['MintFlow_Generated_Xint_plus_Xmic'] = (
            np.asarray(adata_pert.obsm['MintFlow_Generated_Xint'], dtype=np.float64)
            + np.asarray(adata_pert.obsm['MintFLow_Generated_Xmic'], dtype=np.float64)
        )
        adata_pert.obsm['MintFLow_np_MCC'] = np.asarray(result_gen_pert['np_MCC'])
        adata_perturbed_path = mode_dir['generated_adatas'] / f'{safe_name(section_name)}__perturbed.h5ad'
        save_adata(adata_pert, adata_perturbed_path)

        list_cell_id_orig = adata_orig.obs['cell_id'].astype(str).tolist()
        list_cell_id_pert = adata_pert.obs['cell_id'].astype(str).tolist()
        np_mcc_orig = np.asarray(adata_orig.obsm['MintFLow_np_MCC'])
        np_mcc_pert = np.asarray(adata_pert.obsm['MintFLow_np_MCC'])

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
                build_mcc_mask_table(adata_pert, section_name, MODE_CONFIGS[mode_name]['perturbation_label'], changed_pert),
            ],
            ignore_index=True,
        )
        mcc_mask_path = mode_dir['mcc_masks'] / f'{safe_name(section_name)}__mcc_masks.tsv.gz'
        save_table(mcc_mask_df, mcc_mask_path)

        results[section_name] = {
            'adata_original_path': str(shared_generated_dir / f'{safe_name(section_name)}__unperturbed.h5ad'),
            'adata_perturbed_path': str(adata_perturbed_path),
            'mcc_mask_path': str(mcc_mask_path),
            'dict_cellID_to_MCCchanged': dict_cell_id_to_mcc_changed,
        }

        del adata_orig
        del adata_pert
        del result_gen_pert
        del np_mcc_orig
        del np_mcc_pert
        collect_memory()

    with gzip.open(mode_dir['results_pickle'], 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved results pickle -> {mode_dir['results_pickle']}")


def main():
    sc.settings.verbosity = 3
    sc.settings.seed = RANDOM_SEED
    print('GPU Available:', torch.cuda.is_available())
    print(f'Run index: {RUN_INDEX}')
    print(f'Random seed: {RANDOM_SEED}')

    device = torch.device('cpu')
    if torch.cuda.is_available():
        print('GPU Name:', torch.cuda.get_device_name(0))
        device = torch.device('cuda')

    workflow_root = WORKFLOW_ROOT
    workflow_root.mkdir(parents=True, exist_ok=True)
    print(f'Workflow output root: {workflow_root}')

    shared_root = workflow_root / 'shared_unperturbed'
    shared_input_dir = shared_root / 'input_adatas'
    shared_generated_dir = shared_root / 'generated_adatas'
    shared_input_dir.mkdir(parents=True, exist_ok=True)
    shared_generated_dir.mkdir(parents=True, exist_ok=True)

    mode_dirs = {}
    for mode_name in MODE_ORDER:
        mode_root = workflow_root / mode_name
        mode_dirs[mode_name] = {
            'root': mode_root,
            'input_adatas': mode_root / 'input_adatas',
            'generated_adatas': mode_root / 'generated_adatas',
            'mcc_masks': mode_root / 'mcc_masks',
            'results_pickle': mode_root / 'results.pkl.gz',
        }

    model = torch.load(MODEL_CHECKPOINT, map_location='cpu', weights_only=False)
    model.to(device)
    print(f'Loaded checkpoint: {MODEL_CHECKPOINT}')

    data_mintflow, dict_all4_configs = build_mintflow_context(
        PATH_ANNDATA, num_training_epochs=NUM_TRAINING_EPOCHS
    )
    print(f'Built MintFlow context from: {PATH_ANNDATA}')

    ref_level4 = load_tcell_level4_reference()
    target_keys, target_df = load_ifn_macrophage_targets()
    save_table(target_df, workflow_root / 'ifn_target_macrophages.tsv.gz')

    base_adatas = load_base_adatas(ref_level4)
    section_plans, plan_summary_df, target_cells_df = build_section_perturbation_plans(base_adatas, target_keys)
    save_table(plan_summary_df, workflow_root / 'section_target_counts.tsv')
    save_table(target_cells_df, workflow_root / 'section_target_cells.tsv.gz')

    shared_rows = []
    for section_name in TARGET_SECTIONS:
        adata_original = base_adatas[section_name].copy()
        save_adata(adata_original, shared_input_dir / f'{safe_name(section_name)}__unperturbed.h5ad')
        shared_rows.append(
            {
                'section': section_name,
                'n_section_cells': int(adata_original.n_obs),
                'n_target_ifn_macrophages': int(section_plans[section_name]['n_target_cells']),
            }
        )
        del adata_original
    save_table(pd.DataFrame(shared_rows), shared_root / 'section_summary.tsv')

    for mode_name in MODE_ORDER:
        prepare_mode_input_adatas(
            mode_name=mode_name,
            base_adatas=base_adatas,
            section_plans=section_plans,
            mode_dir=mode_dirs[mode_name],
        )

    del base_adatas
    collect_memory()

    generate_shared_unperturbed(
        section_names=TARGET_SECTIONS,
        shared_input_dir=shared_input_dir,
        shared_generated_dir=shared_generated_dir,
        model=model,
        data_mintflow=data_mintflow,
        dict_all4_configs=dict_all4_configs,
        device=device,
    )

    for mode_name in MODE_ORDER:
        print(f"\n{'#' * 80}")
        print(f'Generating mode: {mode_name}')
        print(f"{'#' * 80}")
        generate_mode_outputs(
            section_names=TARGET_SECTIONS,
            mode_name=mode_name,
            mode_dir=mode_dirs[mode_name],
            shared_generated_dir=shared_generated_dir,
            model=model,
            data_mintflow=data_mintflow,
            dict_all4_configs=dict_all4_configs,
            device=device,
        )


if __name__ == '__main__':
    main()
