# CellCharter Node-Removal Reproducibility

This folder reproduces CellCharter neighborhood aggregation without center/self (`include_self=False`).

## Files

- `run_cellcharter_node_removal.py`: main runner.
- `bsub_cellcharter_node_removal.lsf`: LSF job script.
- `requirements.cellcharter-venv.txt`: full package list from CellCharter environment.
- `requirements.cellcharter-only.txt`: minimal CellCharter-focused package list.

## Inputs

Use dataset names with:

- `<RECON_ROOT>/<dataset>/combined.h5ad`

Each input must include `obsm["spatial"]`.

## Outputs

Outputs are written under:

- `<OUT_ROOT>/<dataset>/cellcharter_no_center.h5ad`

## Run locally

```bash
python run_cellcharter_node_removal.py \
  --root /path/to/reconstructed_h5ad \
  --datasets xenium_old_beacon \
  --out-root /path/to/cellcharter_no_center_results
```

With scVI embedding:

```bash
python run_cellcharter_node_removal.py \
  --root /path/to/reconstructed_h5ad \
  --datasets xenium_old_beacon \
  --use-rep X_scVI \
  --sample-key sample
```

## Run with LSF

```bash
bsub < bsub_cellcharter_node_removal.lsf
```

