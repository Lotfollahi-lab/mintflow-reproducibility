# Linear NCEM Reproducibility (Combined Datasets)

This folder contains reproducible scripts to run **linear NCEM** for spatial domain recovery on combined datasets.

## Files

- `run_linear_ncem_combined.py`: main runner.
- `bsub_linear_ncem_combined.lsf`: LSF job template.

## Inputs

Expected input layout under `--in_root`:

- `<in_root>/visiumhd_crc/combined.h5ad`
- `<in_root>/slideseq_kidney/combined.h5ad`
- `<in_root>/xenium_old_kidney/combined.h5ad`
- `<in_root>/xenium_old_beacon/combined.h5ad`
- `<in_root>/xenium_new_psoriasis/combined.h5ad`

Each input must have:

- `obsm["spatial"]`
- dataset-specific cell type key (hardcoded in script)

If `obsp["spatial_connectivities"]` is missing, it is built from `obsm["spatial"]`.

## Outputs

Written to `--out_root/<dataset>/combined_linear_ncem_full.h5ad`.

Main NCEM outputs in each h5ad:

- `obsm["ncem_dmat"]` (design matrix)
- `obsm["ncem_niche"]` (neighborhood composition)

## Run locally

```bash
python run_linear_ncem_combined.py \
  --in_root /lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad \
  --out_root /lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/NCEM/linear_results_parallel
```

Optional: run selected datasets only

```bash
python run_linear_ncem_combined.py \
  --datasets visiumhd_crc,xenium_old_beacon \
  --in_root /path/to/reconstructed_h5ad \
  --out_root /path/to/linear_results_parallel
```

## Run with LSF

```bash
bsub < bsub_linear_ncem_combined.lsf
```

Adjust environment activation and paths in `bsub_linear_ncem_combined.lsf` for your system.
