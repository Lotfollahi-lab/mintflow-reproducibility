# UTAG Main Reproducibility

This folder reproduces UTAG benchmark runs on selected datasets.

## Files

- `run_utag_combined.py`: main UTAG runner.
- `bsub_utag_combined.lsf`: LSF job script.
- `requirements.utag-venv.txt`: full package list from UTAG environment.
- `requirements.utag-only.txt`: minimal UTAG-focused package list.

## Inputs

Use dataset names with:

- `<RECON_ROOT>/<dataset>/combined.h5ad`

Default datasets:

- `xenium_old_kidney` (slide key `section`)
- `xenium_old_beacon` (slide key `sample`)

## Outputs

Default behavior writes `combined_utag.h5ad` next to each input.

You can set `--out-root` to write to:

- `<OUT_ROOT>/<dataset>/combined_utag.h5ad`

## Run locally

```bash
python run_utag_combined.py \
  --root /path/to/reconstructed_h5ad \
  --out-root /path/to/utag_results
```

## Run with LSF

```bash
bsub < bsub_utag_combined.lsf
```

