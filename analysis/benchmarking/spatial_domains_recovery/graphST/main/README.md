# GraphST Main Reproducibility

This folder reproduces GraphST benchmark runs using the existing GraphST runner.

## Files

- `run_graphst_combined.py`: wrapper around the maintained GraphST runner.
- `bsub_graphst_combined.lsf`: LSF job script.
- `requirements.graphst-venv.txt`: full package list from GraphST environment.
- `requirements.graphst-only.txt`: minimal GraphST-focused package list.

## Inputs

Use dataset names with:

- `<RECON_ROOT>/<dataset>/combined.h5ad`

or run a single dataset with `--dataset`.

Each input must include `obsm["spatial"]`.

## Outputs

Outputs are written to:

- `<OUT_ROOT>/<dataset>/<dataset>_graphst.h5ad`
- `<OUT_ROOT>/<dataset>/emb_graphst.npy`

## Run locally

```bash
python run_graphst_combined.py \
  --python /path/to/graphst-env/bin/python \
  --root /path/to/reconstructed_h5ad \
  --out-root /path/to/graphst_results \
  --device cpu \
  --skip-existing
```

## Run with LSF

```bash
bsub < bsub_graphst_combined.lsf
```

Before submitting, edit `bsub_graphst_combined.lsf` for your environment path if needed.

