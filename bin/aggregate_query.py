#!/usr/bin/env python3

"""
Concatenate all sample h5ads for one study, apply gene mapping,
relabeling, label aggregation, and subsample subsample_ref cells per cell type.
"""

import argparse
import os
import random

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

import utils
from utils import map_genes, relabel, aggregate_labels


def subsample_per_celltype(adata, celltype_key, n, seed=42):
    """Sample up to n cells from each cell type."""
    random.seed(seed)
    np.random.seed(seed)
    celltypes = adata.obs[celltype_key].unique()
    keep_idx = []
    for ct in celltypes:
        idx = np.where(adata.obs[celltype_key] == ct)[0]
        if len(idx) > n:
            idx = np.random.choice(idx, size=n, replace=False)
        keep_idx.extend(idx.tolist())
    return adata[keep_idx].copy()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Aggregate sample-level h5ads into a single study-level h5ad."
    )
    parser.add_argument("--study_key", type=str, required=True)
    parser.add_argument("--h5ad_files", type=str, nargs="+", required=True)
    parser.add_argument("--relabel_path", type=str, required=True)
    parser.add_argument("--gene_mapping", type=str, required=True)
    parser.add_argument("--mapping_file", type=str, required=True)
    parser.add_argument(
        "--ref_keys",
        type=str,
        nargs="+",
        default=["subclass", "class", "family", "global"],
    )
    parser.add_argument(
        "--subsample_ref",
        type=int,
        default=500,
        help="Max cells per cell type to keep in the aggregated reference.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--remove_unknown",
        action="store_true",
        help="Remove cells whose most-granular label is 'unknown'.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    random.seed(args.seed)
    np.random.seed(args.seed)
    scvi.settings.seed = args.seed

    gene_mapping = pd.read_csv(args.gene_mapping, sep="\t", header=0)
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")

    # Load and gene-map each sample h5ad
    adatas = []
    sample_names = []
    for f in args.h5ad_files:
        a = ad.read_h5ad(f)
        a = map_genes(a, gene_mapping)
        sample_names.append(os.path.basename(f).replace(".h5ad", ""))
        adatas.append(a)

    # Concatenate across samples
    combined = ad.concat(
        adatas, join="outer", label="sample_id", keys=sample_names
    )
    combined.obs_names_make_unique()

    # Relabel using study-level relabel TSV
    combined = relabel(query=combined, relabel_path=args.relabel_path, sep="\t")

    # Aggregate to higher label levels
    combined.obs = aggregate_labels(
        query=combined.obs,
        mapping_df=mapping_df,
        ref_keys=args.ref_keys,
        predicted=False,
    )

    # Remove unknown labels at most-granular level
    if args.remove_unknown:
        combined = combined[combined.obs[args.ref_keys[0]] != "unknown"].copy()

    # Drop cells where the most-granular label is NaN (not in relabel file)
    combined = combined[combined.obs[args.ref_keys[0]].notna()].copy()

    # Subsample per cell type
    if args.subsample_ref:
        combined = subsample_per_celltype(
            combined, args.ref_keys[0], args.subsample_ref, args.seed
        )

    combined.write_h5ad(f"{args.study_key}.h5ad")
    print(
        f"Wrote {combined.n_obs} cells × {combined.n_vars} genes → {args.study_key}.h5ad"
    )


if __name__ == "__main__":
    main()
