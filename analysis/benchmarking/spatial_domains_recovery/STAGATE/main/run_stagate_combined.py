#!/usr/bin/env python3
import argparse
import inspect
import os

import anndata as ad
import numpy as np


DEFAULT_DATASETS = [
    "slideseq_kidney",
    "visiumhd_crc",
    "xenium_new_psoriasis",
    "xenium_old_beacon",
    "xenium_old_kidney",
]

RAD_CUTOFFS = {
    "slideseq_kidney": 50,
    "visiumhd_crc": 40,
    "xenium_new_psoriasis": 30,
    "xenium_old_beacon": 30,
    "xenium_old_kidney": 30,
}


def import_stagate():
    from STAGATE_pyG import Cal_Spatial_Net, Train_STAGATE

    train_fn = Train_STAGATE
    if not callable(train_fn):
        train_fn = getattr(Train_STAGATE, "train_STAGATE", None)
    if not callable(train_fn):
        raise TypeError("STAGATE_pyG.Train_STAGATE has no callable train function")
    return Cal_Spatial_Net, train_fn


def resolve_dataset(spec: str, root: str):
    if spec.endswith(".h5ad") or os.path.isabs(spec):
        h5ad_path = os.path.abspath(spec)
        base = os.path.basename(h5ad_path)
        dataset = os.path.basename(os.path.dirname(h5ad_path)) if base == "combined.h5ad" else os.path.splitext(base)[0]
        return dataset, h5ad_path
    dataset = spec
    return dataset, os.path.join(root, dataset, "combined.h5ad")


def run_one(dataset, h5ad_path, out_root, rad_cutoff, n_epochs, lr, device, cal_spatial_net, train_stagate):
    print(f"\n=== STAGATE (main): {dataset} ===")
    print("input:", h5ad_path)
    if not os.path.isfile(h5ad_path):
        raise FileNotFoundError(h5ad_path)

    adata = ad.read_h5ad(h5ad_path)
    if "spatial" not in adata.obsm:
        raise ValueError("Missing adata.obsm['spatial'] required for STAGATE.")

    cal_spatial_net(adata, rad_cutoff=rad_cutoff)

    sig = inspect.signature(train_stagate)
    kwargs = {}
    if "n_epochs" in sig.parameters:
        kwargs["n_epochs"] = n_epochs
    if "lr" in sig.parameters:
        kwargs["lr"] = lr
    if "device" in sig.parameters:
        try:
            import torch

            kwargs["device"] = torch.device(device)
        except Exception:
            kwargs["device"] = device

    def run_train(device_override=None):
        run_kwargs = dict(kwargs)
        if device_override is not None and "device" in sig.parameters:
            try:
                import torch

                run_kwargs["device"] = torch.device(device_override)
            except Exception:
                run_kwargs["device"] = device_override
        return train_stagate(adata, **run_kwargs)

    try:
        adata = run_train()
    except RuntimeError as exc:
        if device.startswith("cuda") and "CUDA out of memory" in str(exc):
            print("CUDA OOM detected; retrying on CPU.")
            adata = run_train(device_override="cpu")
        else:
            raise

    if "STAGATE" not in adata.obsm:
        raise RuntimeError(f"{dataset}: STAGATE finished but adata.obsm['STAGATE'] missing")

    out_dir = os.path.join(out_root, dataset)
    os.makedirs(out_dir, exist_ok=True)
    emb = np.asarray(adata.obsm["STAGATE"])
    np.save(os.path.join(out_dir, "emb_stagate.npy"), emb)
    adata.write_h5ad(os.path.join(out_dir, "stagate_full.h5ad"))
    print("saved:", out_dir, "emb shape:", emb.shape)


def main():
    parser = argparse.ArgumentParser(description="Run original STAGATE on selected combined.h5ad datasets.")
    parser.add_argument(
        "--root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/reconstructed_h5ad",
        help="Root containing dataset folders with combined.h5ad.",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=DEFAULT_DATASETS,
        help="Dataset names or explicit .h5ad paths.",
    )
    parser.add_argument(
        "--out-root",
        default="/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/STAGATE/results",
        help="Output root.",
    )
    parser.add_argument("--rad-cutoff", type=int, default=0, help="Global rad cutoff override. 0 = per-dataset defaults.")
    parser.add_argument("--n-epochs", type=int, default=1000, help="Training epochs.")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate.")
    parser.add_argument("--device", default="cuda", help="Device (cuda/cpu).")
    args = parser.parse_args()

    os.makedirs(args.out_root, exist_ok=True)
    cal_spatial_net, train_stagate = import_stagate()

    for spec in args.datasets:
        dataset, h5ad_path = resolve_dataset(spec, os.path.abspath(args.root))
        rad_cutoff = args.rad_cutoff if args.rad_cutoff > 0 else RAD_CUTOFFS.get(dataset, 50)
        run_one(
            dataset,
            h5ad_path,
            args.out_root,
            rad_cutoff,
            args.n_epochs,
            args.lr,
            args.device,
            cal_spatial_net,
            train_stagate,
        )


if __name__ == "__main__":
    main()

