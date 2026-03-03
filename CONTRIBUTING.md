# Contributing to MintFlow Reproducibility

Thank you for contributing your analysis to the MintFlow reproducibility repository. This guide helps co-authors add their analyses in a consistent, reviewer-friendly format.

## How to Add Your Analysis

### 1. Create an Analysis Subfolder

Create a new folder under `analysis/` with a descriptive name:

- **By figure**: `analysis/figure5_application_name/`
- **By application**: `analysis/breast_cancer/`, `analysis/lung_atlas/`
- Use lowercase with underscores (e.g., `figure2_3_eczema`)

### 2. Add Your Files

Place your notebooks (`.ipynb`), scripts (`.py`, `.R`), and any small helper files in the subfolder.

**Naming convention**: Use clear, descriptive names such as:
- `01_data_prep.ipynb`
- `02_main_analysis.ipynb`
- `Fig7_survival_analysis.R`

### 3. Add a README

Create a `README.md` in your subfolder with:

- **Brief description** of the analysis and which figure(s) it produces
- **Table of files** with short descriptions
- **Data requirements**: paths, formats, or links to external data
- **Dependencies**: any packages beyond the main environment (e.g., R packages)

Example:

```markdown
# Figure X: Application Name

Brief description of what this analysis does.

| File | Description |
|------|-------------|
| `01_prep.ipynb` | Data preprocessing |
| `02_analysis.ipynb` | Main analysis producing Figure X |

## Data

- Input: `path/to/adata.h5ad` (or GDrive link)
- Output: Figures Xa, Xb
```

### 4. Use Portable Paths

- Prefer **relative paths** where possible
- Document **absolute paths** or environment variables if data lives elsewhere
- Add a cell or comment at the top explaining how to set `DATA_DIR` or similar

### 5. Environment

- If your analysis needs **extra packages**, add them to `envs/environment.yaml` or document them in your README
- For **R-only** analyses, list required R packages in your README

### 6. Submit

1. Commit your changes
2. Push to your branch
3. Open a Pull Request, or coordinate with the maintainers for direct merge

## Checklist for New Analyses

- [ ] New subfolder under `analysis/`
- [ ] `README.md` with description and file table
- [ ] Data requirements documented
- [ ] Paths are documented or configurable
- [ ] Notebooks run (or known issues noted)

## Questions?

Open an issue or contact the maintainers.
