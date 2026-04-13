# Reproducibility (data and model weights)

Reviewers expect each analysis to state **where to obtain inputs** (data) and **trained model checkpoints** (weights), not only the code. This repository keeps large files off Git; instead we document **stable download links** and **versions** next to each analysis.

## What each analysis must document

Every folder under `analysis/<name>/` should include a `README.md` with a **Reproducibility** block (copy the template below). Optionally add the same block as a short markdown cell at the top of each notebook.

### 1. Data

| Item | What to record |
|------|----------------|
| **Source** | Zenodo DOI, Figshare, institutional download, or GDrive folder URL |
| **Files** | Exact filenames or patterns (e.g. `adata_eczema.h5ad`) |
| **Version** | Snapshot date, record ID, or checksum (SHA256) if available |

### 2. Model weights (MintFlow / other models)

| Item | What to record |
|------|----------------|
| **Checkpoint** | URL or Zenodo file for the `.pt` / `.ckpt` / `saved model` used in the paper |
| **Training config** | Commit hash of training code, or config filename bundled with the checkpoint |
| **MintFlow version** | `pip show mintflow` or git tag of [mintflow](https://github.com/Lotfollahi-lab/mintflow) used to train or load weights |

If weights are not public yet, state **“available on request”** or the **planned Zenodo record** and update when live.

### 3. Software environment

- Conda: `envs/environment.yaml` (pin versions where critical).
- One-line note in the analysis README if a notebook needs a **different** env file.

## Where to host large files

| Option | Best for |
|--------|----------|
| **[Zenodo](https://zenodo.org)** | DOI, versioning, reviewer-friendly citation |
| **Figshare / institutional repository** | Same idea as Zenodo |
| **Google Drive** | OK if you add a **stable folder link** and file list in README |
| **Git LFS** | Only for smaller weights; avoid multi-GB files in Git |

### Restricted or password-protected data

**Zenodo is for public releases** (anyone with the DOI can download). If you must **not** share data openly, do not rely on Zenodo as the primary host for those files.

Use instead:

- **Peer review:** A private folder (Google Drive, Dropbox, Box, institutional storage) with **password or invite-only access**; share credentials with the journal or reviewers as **confidential supplementary material** (check journal policy).
- **After publication:** State **“available to qualified researchers on reasonable request”** (or access via **dbGaP / EGA** for human genomics), and document the **request process** (email, committee) in the manuscript and in each analysis `README.md`.
- **Optional later:** Deposit a **public subset** or **metadata-only** record on Zenodo (DOI for “what exists”) while full data stay restricted elsewhere.

In analysis READMEs, write e.g. **“Data: restricted — contact [lab email] / access via [process]. Not on public Zenodo.”** so expectations are clear for reviewers and readers.

## Central index (optional)

Maintain [`datasets/README.md`](datasets/README.md) as a **single table** listing all manuscript datasets and model bundles with DOIs/links, so reviewers see everything in one place.

## Template (paste into `analysis/<your_folder>/README.md`)

```markdown
## Reproducibility

### Data
- **Download:** [Zenodo DOI or URL](https://...)
- **Files used:** `file1.h5ad`, `file2.csv`
- **Checksum (optional):** `sha256:...`

### Model weights
- **MintFlow (or other) checkpoint:** [URL or Zenodo](https://...)
- **Compatible package version:** `mintflow==x.y.z` (or git commit `abc123`)

### Environment
- `conda env create -f ../../envs/environment.yaml`
```

## Notebooks

At the top of each notebook, use a small **configuration cell** so paths are obvious:

```python
# Reproducibility — set before running
DATA_DIR = "/path/to/downloaded_data"  # see analysis/<folder>/README.md
MODEL_CKPT = "/path/to/mintflow_weights.pt"  # same README
```

Replace with `os.environ.get("MINTFLOW_DATA_DIR", "...")` if you prefer environment variables.
