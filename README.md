# MintFlow Reproducibility

This repository contains the code to reproduce the analyses and benchmarking experiments performed in the MintFlow [manuscript](https://www.biorxiv.org/content/10.1101/2025.06.24.661094v1). The MintFlow source code can be found [here](https://github.com/Lotfollahi-lab/mintflow).

This structure is designed for reviewer readability and to enable co-authors to contribute their analyses. See [CONTRIBUTING.md](CONTRIBUTING.md) for how to add your analysis.

## Repository Structure

```
mintflow-reproducibility/
├── analysis/           # Reproducible analyses organized by figure/application
│   ├── figure2_3_eczema/    # Eczema scRNA-seq and drug2cell analyses
│   ├── figure4_melanoma/    # Melanoma spatial transcriptomics
│   └── kidney_cancer/       # Renal cell carcinoma (RCC) analyses
├── datasets/           # Dataset documentation and download info
├── envs/               # Conda environment specifications
├── utils/              # Shared utility functions
└── README.md
```

## Installation

### Standard

1. Clone the mintflow-reproducibility repository and navigate into it:
   ```bash
   git clone https://github.com/Lotfollahi-lab/mintflow-reproducibility.git
   cd mintflow-reproducibility
   ```

2. (Optional) Install the Libmamba solver to make the installation faster:
   ```bash
   conda update -n base conda
   conda install -n base conda-libmamba-solver
   conda config --set solver libmamba
   ```

3. Create the mintflow-reproducibility conda environment:
   ```bash
   conda env create -f envs/environment.yaml
   conda activate mintflow-reproducibility
   ```

4. Install MintFlow (if not already installed):
   ```bash
   pip install mintflow
   ```
   Or install from source: `pip install git+https://github.com/Lotfollahi-lab/mintflow.git`

### R Analyses

The kidney cancer survival analysis (`Fig6&S14_MintFlow_RCCanalysis_TGCASurvival.R`) requires R. Install R and the required packages as indicated in the script.

## Data & Models

Preprocessed data and **trained model weights** are linked from each analysis folder and summarized in [`datasets/README.md`](datasets/README.md). Large assets are hosted externally (e.g. Zenodo with a DOI, or GDrive) — see [**Reproducibility**](REPRODUCIBILITY.md) for what reviewers expect and a README template.

When running analyses, set paths as documented in each `analysis/<name>/README.md` (or the configuration cell at the top of each notebook).

## Contributing Your Analysis

We welcome contributions from co-authors. **Where to add your analysis:**

| Add your analysis here | Example |
|------------------------|---------|
| `analysis/<your_folder>/` | `analysis/figure5_application/` or `analysis/breast_cancer/` | 
| Add a `README.md` inside your folder | Describes files and data requirements |

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full guide (naming, paths, checklist).

## Reference

```
@article{Akbarnejad2025,
  author    = {Akbarnejad, A. et al.},
  title     = {Mapping and reprogramming microenvironment-induced cell states in human disease using generative AI},
  journal   = {bioRxiv},
  year      = {2025},
  doi       = {10.1101/2025.06.24.661094},
  url       = {https://doi.org/10.1101/2025.06.24.661094}
}
```
