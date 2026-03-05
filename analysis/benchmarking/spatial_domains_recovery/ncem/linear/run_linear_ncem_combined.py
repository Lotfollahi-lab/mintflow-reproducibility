#!/usr/bin/env python3
import argparse
import os
import re
import sys

import anndata as ad
import ncem.tl as tl
import numpy as np
import scipy.sparse as sp

sys.setrecursionlimit(200000)

KEY_GRAPH = "spatial_connectivities"
SEED = 0
N_NEIGHBORS_DEFAULT = 10

DATASETS = {
    "visiumhd_crc": {
        "key_type": "DeconvolutionLabel2",
        "min_cells_per_type": 50,
        "n_neighbors": 5,
    },
    "slideseq_kidney": {
        "key_type": "scANVI_pred",
        "min_cells_per_type": 20,
        "n_neighbors": N_NEIGHBORS_DEFAULT,
    },
    "xenium_old_kidney": {
        "key_type": "level_2_cell_type",
        "min_cells_per_type": 20,
        "n_neighbors": N_NEIGHBORS_DEFAULT,
    },
    "xenium_old_beacon": {
        "key_type": "broad_celltypes",
        "min_cells_per_type": 20,
        "n_neighbors": N_NEIGHBORS_DEFAULT,
    },
    "xenium_new_psoriasis": {
        "key_type": "lvl2_annotation",
        "min_cells_per_type": 2,
        "n_neighbors": N_NEIGHBORS_DEFAULT,
    },
}


def safe_token(s: str) -> str:
    s = str(s)
    s = re.sub(r"[^0-9a-zA-Z_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    if not s:
        s = "NA"
    if s[0].isdigit():
        s = "ct_" + s
    return s


def patch_disable_wald_tests() -> None:
    import ncem.tl.fit.backend.linear_model as lm

    def _no_test_standard(*args, **kwargs):
        return kwargs.get("adata", args[0] if args else None)

    lm.test_standard = _no_test_standard
    print("[patch] lm.test_standard -> no-op")


def ensure_spatial_graph(adata: ad.AnnData, n_neighbors: int) -> None:
    if KEY_GRAPH in adata.obsp:
        return
    if "spatial" not in adata.obsm:
        raise RuntimeError("Missing obsm['spatial']; cannot build spatial graph")

    import scanpy as sc

    print(f"[graph] building {KEY_GRAPH} with n_neighbors={n_neighbors}")
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep="spatial",
        key_added="spatial",
        random_state=SEED,
    )
    if KEY_GRAPH not in adata.obsp:
        raise RuntimeError(f"Failed to build {KEY_GRAPH}")


def run_one(dataset: str, in_root: str, out_root: str) -> None:
    cfg = DATASETS[dataset]
    in_path = os.path.join(in_root, dataset, "combined.h5ad")
    if not os.path.exists(in_path):
        print(f"[skip] missing {in_path}")
        return

    out_dir = os.path.join(out_root, dataset)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "combined_linear_ncem_full.h5ad")

    print(f"\\n=== {dataset} / combined ===")
    print("load:", in_path)
    adata = ad.read_h5ad(in_path)
    print("shape:", adata.shape)

    key_type_raw = cfg["key_type"]
    key_type_san = key_type_raw + "_sanitized"
    adata.obs[key_type_san] = adata.obs[key_type_raw].map(safe_token).astype("category")
    adata.obs[key_type_san] = adata.obs[key_type_san].cat.remove_unused_categories()

    vc = adata.obs[key_type_san].value_counts()
    keep = vc[vc >= int(cfg["min_cells_per_type"])].index
    adata = adata[adata.obs[key_type_san].isin(keep)].copy()
    adata.obs[key_type_san] = adata.obs[key_type_san].astype("category").cat.remove_unused_categories()
    print("after type filter:", adata.shape, "n_types:", adata.obs[key_type_san].nunique())

    ensure_spatial_graph(adata, n_neighbors=int(cfg.get("n_neighbors", N_NEIGHBORS_DEFAULT)))
    if not sp.issparse(adata.obsp[KEY_GRAPH]):
        raise RuntimeError(f"{KEY_GRAPH} is not sparse")

    patch_disable_wald_tests()
    tl.linear_ncem(
        adata=adata,
        key_type=key_type_san,
        key_graph=KEY_GRAPH,
        formula="~0",
    )
    print("[ok] linear_ncem done")

    adata.write_h5ad(out_path)
    print("saved:", out_path)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--in_root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad",
    )
    ap.add_argument(
        "--out_root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/NCEM/linear_results_parallel",
    )
    ap.add_argument(
        "--datasets",
        default=",".join(DATASETS.keys()),
        help="Comma-separated dataset names",
    )
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_root, exist_ok=True)
    datasets = [x.strip() for x in args.datasets.split(",") if x.strip()]

    unknown = [d for d in datasets if d not in DATASETS]
    if unknown:
        raise ValueError(f"Unknown datasets: {unknown}")

    for ds in datasets:
        run_one(ds, args.in_root, args.out_root)


if __name__ == "__main__":
    main()
