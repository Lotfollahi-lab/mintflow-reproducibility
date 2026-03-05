#!/usr/bin/env python3
import argparse
import os

import scanpy as sc
import squidpy as sq
import cellcharter as cc


DEFAULT_DATASETS = ["xenium_old_beacon"]


def main():
    parser = argparse.ArgumentParser(description="Run CellCharter neighborhood aggregation without center/self.")
    parser.add_argument(
        "--root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad",
        help="Root directory containing <dataset>/combined.h5ad.",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=DEFAULT_DATASETS,
        help="Dataset names under --root.",
    )
    parser.add_argument(
        "--out-root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/cellcharter/results_no_center",
        help="Output root directory.",
    )
    parser.add_argument("--n-layers", type=int, default=3, help="CellCharter n_layers.")
    parser.add_argument("--use-rep", default="", help="Optional obsm key, e.g. X_scVI.")
    parser.add_argument("--sample-key", default="", help="Optional sample key in obs.")
    parser.add_argument("--out-key", default="X_cellcharter", help="Output key in obsm.")
    args = parser.parse_args()

    os.makedirs(args.out_root, exist_ok=True)
    for ds in args.datasets:
        in_path = os.path.join(args.root, ds, "combined.h5ad")
        if not os.path.exists(in_path):
            raise FileNotFoundError(in_path)

        print(f"\n=== CellCharter node-removal: {ds} ===")
        adata = sc.read_h5ad(in_path)
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

        kwargs = dict(n_layers=args.n_layers, include_self=False, out_key=args.out_key)
        if args.use_rep:
            kwargs["use_rep"] = args.use_rep
        if args.sample_key:
            kwargs["sample_key"] = args.sample_key
        cc.gr.aggregate_neighbors(adata, **kwargs)

        out_dir = os.path.join(args.out_root, ds)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "cellcharter_no_center.h5ad")
        adata.write_h5ad(out_path)
        print("saved:", out_path)


if __name__ == "__main__":
    main()

