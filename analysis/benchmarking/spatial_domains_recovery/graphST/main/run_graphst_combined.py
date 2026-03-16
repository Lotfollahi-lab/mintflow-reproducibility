#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Wrapper to run GraphST combined benchmark.")
    parser.add_argument(
        "--python",
        default="/software/cellgen/team361/ms83/envs/graphst-env/bin/python",
        help="Python interpreter for GraphST environment.",
    )
    parser.add_argument(
        "--runner",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/GraphST/run_graphst_all_combined.py",
        help="Path to the GraphST runner.",
    )
    parser.add_argument(
        "--root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad",
        help="Root with <dataset>/combined.h5ad inputs.",
    )
    parser.add_argument(
        "--out-root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/GraphST/results",
        help="Output root for GraphST results.",
    )
    parser.add_argument(
        "--log-dir",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/GraphST/logs",
        help="Per-dataset log directory used by the underlying runner.",
    )
    parser.add_argument("--device", default="cpu", help="Device for GraphST runner.")
    parser.add_argument("--n-neighbors", type=int, default=10, help="KNN neighbors for sparse graph mode.")
    parser.add_argument("--dataset", default=None, help="Optional single dataset name.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip datasets with existing outputs.")
    args = parser.parse_args()

    cmd = [
        args.python,
        args.runner,
        "--root",
        args.root,
        "--out-root",
        args.out_root,
        "--log-dir",
        args.log_dir,
        "--device",
        args.device,
        "--sparse-graph",
        "--n-neighbors",
        str(args.n_neighbors),
    ]
    if args.dataset:
        cmd += ["--dataset", args.dataset]
    if args.skip_existing:
        cmd += ["--skip-existing"]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

