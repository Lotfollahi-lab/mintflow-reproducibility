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

Create a local environment in this folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.stagate-only.txt
# if needed for exact reproducibility, use:
# pip install -r requirements.stagate-venv.txt
```

## Inputs

Use either:

- explicit `.h5ad` paths via `--datasets`, or
- dataset names with `--root` where each dataset has `combined.h5ad`.

Suggested layout for named datasets:

- `<RECON_ROOT>/<dataset>/combined.h5ad`

Each input must include `obsm["spatial"]`.

## Outputs

By default, outputs are written under:

- `<OUTPUT_ROOT>/<dataset>/`

Per dataset:

- `emb_stagate.npy`
- `stagate_full.h5ad`

## Run locally

```bash
source .venv/bin/activate
python run_stagate_combined_node_removal.py
```

Optional override:

```bash
python run_stagate_combined_node_removal.py \
  --datasets /path/to/sample1.h5ad /path/to/sample2.h5ad \
  --out-root /path/to/stagate_node_removal_results \
  --device cuda
```

## Run with LSF

```bash
bsub < bsub_stagate_node_removal.lsf
```

Before submitting, edit `bsub_stagate_node_removal.lsf` and set `PY` to your environment's Python path.
