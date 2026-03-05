# STAGATE Main Reproducibility

This folder reproduces the **original STAGATE** benchmark runs (center/self contribution intact).

## Files

- `run_stagate_combined.py`: main runner.
- `bsub_stagate_combined.lsf`: LSF job script.
- `requirements.stagate-venv.txt`: full `pip freeze` from the STAGATE environment used.
- `requirements.stagate-only.txt`: minimal STAGATE-focused package list.

## Environment

Runner uses:

- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/STAGATE/.venv_stagate/bin/python`

## Inputs

Default root:

- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad`

Default datasets:

- `slideseq_kidney`
- `visiumhd_crc`
- `xenium_new_psoriasis`
- `xenium_old_beacon`
- `xenium_old_kidney`

Each dataset expects `<dataset>/combined.h5ad` with `obsm["spatial"]`.

## Outputs

Default output root:

- `/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/STAGATE/results`

Per dataset:

- `emb_stagate.npy`
- `stagate_full.h5ad`

## Run locally

```bash
python run_stagate_combined.py
```

Optional:

```bash
python run_stagate_combined.py --datasets visiumhd_crc xenium_old_beacon --device cuda --rad-cutoff 40
```

## Run with LSF

```bash
bsub < bsub_stagate_combined.lsf
```

