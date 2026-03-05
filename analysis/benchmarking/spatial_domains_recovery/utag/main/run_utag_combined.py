#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import scanpy as sc
from utag import utag


DEFAULT_DATASETS = [
    ("xenium_old_kidney", "section"),
    ("xenium_old_beacon", "sample"),
]


def main():
    parser = argparse.ArgumentParser(description="Run UTAG on selected combined.h5ad datasets.")
    parser.add_argument(
        "--root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad",
        help="Root directory containing <dataset>/combined.h5ad.",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=[x[0] for x in DEFAULT_DATASETS],
        help="Dataset names to run.",
    )
    parser.add_argument(
        "--max-dist",
        type=float,
        default=20.0,
        help="UTAG max_dist parameter.",
    )
    parser.add_argument(
        "--apply-umap",
        action="store_true",
        help="Enable UMAP in UTAG.",
    )
    parser.add_argument(
        "--out-root",
        default="",
        help="Optional output root; if empty, writes next to input combined.h5ad.",
    )
    args = parser.parse_args()

    slide_key_map = dict(DEFAULT_DATASETS)

    for ds in args.datasets:
        if ds not in slide_key_map:
            raise ValueError(f"{ds}: missing slide_key mapping; add it to DEFAULT_DATASETS.")

        in_path = Path(args.root) / ds / "combined.h5ad"
        if not in_path.exists():
            raise FileNotFoundError(in_path)

        print(f"\n=== UTAG: {ds} ===")
        print("input:", in_path)
        adata = sc.read_h5ad(in_path)
        adata = utag(adata, slide_key=slide_key_map[ds], max_dist=args.max_dist, apply_umap=args.apply_umap)

        if args.out_root:
            out_dir = Path(args.out_root) / ds
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / "combined_utag.h5ad"
        else:
            out_path = in_path.with_name("combined_utag.h5ad")

        adata.write_h5ad(out_path)
        print("saved:", out_path)


if __name__ == "__main__":
    main()

