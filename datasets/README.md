# Datasets and model weights

Large files (`.h5ad` data, `.pkl` model weights) stay **outside Git**. They are hosted in a private Google Drive folder (password shared with reviewers via the journal).

## Registry

The **single source of truth** for all data and weights is:

- **[`data_and_weights_registry.xlsx`](data_and_weights_registry.xlsx)** (Excel file in this folder)

Each row maps a notebook to the exact files it needs (data or weights), with download links.

### How to use

1. Open the Excel file
2. Find the row(s) for the notebook you want to run
3. Download the listed files from the GDrive link
4. Place them in the path shown (or set `DATA_DIR` / `MODEL_CKPT` at the top of the notebook)

### For co-authors adding files

Add one row per file (`.h5ad` or `.pkl`) to the Excel sheet:

| Column | What to fill in |
|--------|----------------|
| `analysis_folder` | Your folder name under `analysis/` |
| `notebook_or_script` | Which notebook needs this file |
| `filename` | Exact filename (e.g. `adata_eczema.h5ad`) |
| `file_type` | `data` or `weights` |
| `description` | Short description |
| `size_approx` | e.g. `2.3 GB` |
| `gdrive_path_or_link` | Subfolder path in GDrive or direct link |
| `notes` | Any extra info (version, access restrictions) |

**Guidance:** See [REPRODUCIBILITY.md](../REPRODUCIBILITY.md) for hosting options and access policies.
