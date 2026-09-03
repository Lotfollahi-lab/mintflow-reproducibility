#!/usr/bin/env python

import gc
import os
import pickle

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import torch
import mintflow
from anndata.io import read_elem

sc.settings.seed = 0
np.random.seed(0)

fname_checkpoint_model = \
'/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/ICB_TCell_Training_unseen_scRNAseq_subtype/outputs/model_epoch_6.pt'
fname_dictall4configs = \
'/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/ICB_TCell_Training_unseen_scRNAseq_subtype/outputs/dict_all4_configs.pkl'
dir_save = \
'/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Final/Perturbations/Akbar_notebook_final/orig_batch09_50seed'
root_generations = f'{dir_save}/generations_50seed'
root_crops = '/lustre/scratch126/cellgen/lotfollahi/zh4/Mintflow/Subclustering_analysis/adata_TLS_cropped_combined'
macrophage_reference_path = '/nfs/team361/cl35/Macs/adata_macrophages.h5ad'

relabel_map = {
    'CD4+ T lymphocyte': 'CD4+ T lymphocyte ICB',
    'CD8+ T lymphocyte': 'CD8+ T lymphocyte ICB',
    'DN T lymphocyte': 'DN T lymphocyte ICB',
}
ct_filter = 'Tumour associated macrophage'
dose_fraction = {'dose_response_25pct': 0.25, 'dose_response_50pct': 0.5,
                 'dose_response_75pct': 0.75}

sections = [
    'CV1-KID-0-FO-1',
    'CV1-KID-0-FT-2',
    'CV5-KID-0-FO-1',
    'CV6-KID-0-FT-1',
    'CV7-KID-0-FT-2-s3',
    'CV9-KID-0-FT-2',
]

modes_to_generate = [
    'unperturbed',
    'full_replacement',
    'dose_response_25pct',
    'dose_response_50pct',
    'dose_response_75pct',
    'random_cells_uniform',
]

batch_idx_generation = 9  # Batch ID of section with highest wasserstein distance
num_realisations_innerloop = 10
num_seeds_outerloop = 50
seeds_per_job = 10


def arm_path(kind, section, mode, seed):
    sub = 'shared_unperturbed' if mode == 'unperturbed' else mode
    tag = 'unperturbed' if mode == 'unperturbed' else 'perturbed'
    return f'{root_generations}/{sub}/{kind}_seed{seed:02d}/{section}__{tag}.h5ad'


def write_h5ad(adata, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path)


### Relabelling helper
def relabel(adata, mode, rng):
    old = adata.obs['level_3_cell_type'].astype(str)
    new = old.copy()
    target = old.isin(list(relabel_map))

    if mode == 'full_replacement':
        new.loc[target] = new.loc[target].map(relabel_map)
    elif mode in dose_fraction:
        for source, dest in relabel_map.items():
            idx = rng.permutation(np.flatnonzero(old.to_numpy() == source))
            new.iloc[idx[:int(np.rint(dose_fraction[mode] * idx.size))]] = dest
    elif mode == 'random_cells_uniform':
        n = int(target.sum())
        donor = old != ct_filter
        selected = rng.choice(np.flatnonzero(donor.to_numpy()), size=n, replace=False)
        new.iloc[selected] = rng.choice(old[donor].unique(), size=n, replace=True)
    elif mode != 'unperturbed':
        raise ValueError(mode)

    adata.obs['level_3_cell_type'] = new.to_numpy()
    return adata


### Tam level4 subtype reference, used by add_level4() below
with h5py.File(macrophage_reference_path, 'r') as f:
    ref_level4 = pd.Series(
        np.asarray(read_elem(f['obs']['Annotation_TAMs_merged'])).astype(str),
        index=pd.MultiIndex.from_arrays(
            [np.asarray(read_elem(f['obs']['batch'])).astype(str),
             np.asarray(read_elem(f['obs']['_index'])).astype(str)],
            names=['batch', 'cell_id'],
        ),
    )
ref_level4 = ref_level4[~ref_level4.index.duplicated(keep='first')]


def add_level4(adata):
    key = pd.MultiIndex.from_arrays(
        [adata.obs['batch'].astype(str).to_numpy(), adata.obs_names.astype(str)],
        names=['batch', 'cell_id'],
    )
    adata.obs['level_4_cell_type'] = pd.Categorical(ref_level4.reindex(key).to_numpy())
    return adata


# this job covers seeds_per_job consecutive seeds, the batch it starts at comes from the jobscript
seed_offset = int(os.environ['SEED_BATCH_OFFSET'])
seeds = list(range(seed_offset, seed_offset + seeds_per_job))
if seeds[0] < 0 or seeds[-1] >= num_seeds_outerloop:
    raise ValueError(
        f'seeds must be in [0, {num_seeds_outerloop}), got {seeds[0]}-{seeds[-1]} '
        f'(SEED_BATCH_OFFSET={seed_offset})'
    )

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f'seeds {seeds[0]}-{seeds[-1]}, device {device}', flush=True)

todo = [(section, mode, seed) for seed in seeds for mode in modes_to_generate
        for section in sections
        if not os.path.exists(arm_path('seed_adatas', section, mode, seed))]
print(f'{len(todo)} generations missing for seeds {seeds[0]}-{seeds[-1]}', flush=True)

if not todo:
    print('nothing to do', flush=True)
    raise SystemExit

checkpoint_mintflow = torch.load(fname_checkpoint_model, map_location='cpu', weights_only=False)
checkpoint_mintflow.to(device)

with open(fname_dictall4configs, 'rb') as f:
    dict_all4_configs = pickle.load(f)
data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
kwargs_neighbourhood_graph = {
    k: {'True': True, 'False': False}.get(v, v)
    for k, v in dict_all4_configs['config_data_train'][0]['config_neighbourhood_graph'].items()
}
print('loaded model and data', flush=True)

base_cache = {}


def load_base(section):
    """The crop, its level4 labels and its graph depend only on the section, so build each once."""
    if section not in base_cache:
        base = add_level4(sc.read_h5ad(f'{root_crops}/{section}_crop.h5ad'))
        base.uns = {}
        if 'spatial_connectivities' in base.obsp:
            del base.obsp['spatial_connectivities']
        sq.gr.spatial_neighbors(adata=base, **kwargs_neighbourhood_graph)
        base_cache[section] = base
    return base_cache[section]


for section, mode, seed in todo:
    torch.manual_seed(seed)
    adata = relabel(load_base(section).copy(), mode, np.random.default_rng(seed))
    write_h5ad(adata, arm_path('label_adatas', section, mode, seed))

    result_generation = mintflow.generate_insilico_ST_data(
        adata=adata,
        obskey_celltype='level_3_cell_type',
        obspkey_neighbourhood_graph='spatial_connectivities',
        device=device,
        batch_index_trainingdata=batch_idx_generation,
        num_generated_realisations=num_realisations_innerloop,
        model=checkpoint_mintflow,
        data_mintflow=data_mintflow,
        dict_all4_configs=dict_all4_configs,
        estimate_spatial_sizefactors_on_sections=[batch_idx_generation],
    )
    realisations = result_generation['list_generated_realisations_ie_expressions']
    xmic = np.stack([r['MintFLow_Generated_Xmic'] for r in realisations]).mean(0)
    xint = np.stack([r['MintFlow_Generated_Xint'] for r in realisations]).mean(0)
    np_mcc = result_generation['np_MCC'] + 0.0

    del result_generation, realisations
    gc.collect()

    adata.obsm['MintFLow_Generated_Xmic'] = xmic
    adata.obsm['MintFlow_Generated_Xint'] = xint
    adata.obsm['MintFlow_Generated_Xint_plus_Xmic'] = xint + xmic
    adata.obsm['MintFLow_np_MCC'] = np_mcc
    write_h5ad(adata, arm_path('seed_adatas', section, mode, seed))
    print(section, mode, seed, adata.n_obs, flush=True)

    del adata, xmic, xint
    gc.collect()
