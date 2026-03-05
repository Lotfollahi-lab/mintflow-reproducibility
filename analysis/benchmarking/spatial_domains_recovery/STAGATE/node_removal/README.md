# STAGATE Node-Removal Reproducibility

This folder reproduces the **node-removal STAGATE variant** where center/self aggregation is removed.

## Files

- `run_stagate_combined_node_removal.py`: main runner.
- `bsub_stagate_node_removal.lsf`: LSF job script.
- `local_stagate/`: patched local STAGATE package copy used for no-self aggregation.
- `requirements.stagate-venv.txt`: full `pip freeze` from the STAGATE environment used.
- `requirements.stagate-only.txt`: minimal STAGATE-focused package list.

## Node-removal detail

The patched code in `local_stagate/STAGATE_pyG/utils.py` removes identity/self-loop addition during PyG graph transfer, so central nodes are not explicitly included in neighborhood aggregation.

## Environment

Runner uses:

- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/STAGATE/.venv_stagate/bin/python`

## Inputs

Default datasets are explicit paths:

- `/nfs/team361/aa36/OnGit/inflow-reproducibility/Analysis/18_Melanoma_Jan20Runs/NonGit/Data/preprocBy_nb1_preprocessing_V2.h5ad`
- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad/xenium_old_beacon/combined.h5ad`
- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad/xenium_old_kidney/combined.h5ad`

Each input must include `obsm["spatial"]`.

## Outputs

Default output root:

- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/STAGATE/without_centre/results`

Per dataset:

- `emb_stagate.npy`
- `stagate_full.h5ad`

## Run locally

```bash
python run_stagate_combined_node_removal.py
```

Optional override:

```bash
python run_stagate_combined_node_removal.py --datasets /path/to/combined.h5ad --device cuda
```

## Run with LSF

```bash
bsub < bsub_stagate_node_removal.lsf
```

