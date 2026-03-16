# STAGATE Main Reproducibility

This folder reproduces the **original STAGATE** benchmark runs (center/self contribution intact).

## Files

- `run_stagate_combined.py`: main runner.
- `bsub_stagate_combined.lsf`: LSF job script.
- `requirements.stagate-venv.txt`: full `pip freeze` from the STAGATE environment used.
- `requirements.stagate-only.txt`: minimal STAGATE-focused package list.

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

Use dataset names with `--root` and this layout:

- `<RECON_ROOT>/<dataset>/combined.h5ad`

or pass explicit `.h5ad` paths with `--datasets`.

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
python run_stagate_combined.py
```

Optional:

```bash
python run_stagate_combined.py \
  --datasets visiumhd_crc xenium_old_beacon \
  --root /path/to/reconstructed_h5ad \
  --out-root /path/to/stagate_results \
  --device cuda \
  --rad-cutoff 40
```

## Run with LSF

```bash
bsub < bsub_stagate_combined.lsf
```

Before submitting, edit `bsub_stagate_combined.lsf` and set `PY` to your environment's Python path.
