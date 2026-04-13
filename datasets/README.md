# Datasets and model weights

Large files (`.h5ad` data, `.pkl` model weights) stay **outside Git**. They are hosted in a private Google Drive folder (password shared with reviewers via the journal).

## Registry

The list of **all data files and model weights** required to reproduce each notebook is maintained in the following Google Sheet:

**[Data & Weights Registry (Google Sheet)](https://docs.google.com/spreadsheets/d/1MZ575oXGSmGBPvi-qBCIgjDxeWn2tQmCxVbO9LWlhNs/edit?usp=sharing)**

Each row maps a notebook to the exact `.h5ad` data file and `.pkl` model weights it needs to run.

### How to use

1. Open the [Google Sheet](https://docs.google.com/spreadsheets/d/1MZ575oXGSmGBPvi-qBCIgjDxeWn2tQmCxVbO9LWlhNs/edit?usp=sharing)
2. Find the row(s) for the notebook you want to run
3. Download the listed data and weight files from the provided links
4. Place them in the path shown (or set `DATA_DIR` / `MODEL_CKPT` at the top of the notebook)

### For co-authors

Edit the Google Sheet directly to add or update entries for your analysis notebooks.

**Guidance:** See [REPRODUCIBILITY.md](../REPRODUCIBILITY.md) for hosting options and access policies.
