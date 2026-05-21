#!/usr/bin/env python3

"""
Take an already-relabeled, aggregated study h5ad and project it
through the Census SCVI model into the shared latent space.
"""

import argparse
import random

import anndata as ad
import numpy as np
import scvi

from utils import process_query


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="SCVI projection for an aggregated study h5ad."
    )
    parser.add_argument("--h5ad_path", type=str, required=True)
    parser.add_argument("--study_key", type=str, required=True)
    parser.add_argument("--model_path", type=str, required=True)
    parser.add_argument("--batch_key", type=str, default="sample_id")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_arguments()

    random.seed(args.seed)
    np.random.seed(args.seed)
    scvi.settings.seed = args.seed

    query = ad.read_h5ad(args.h5ad_path)

    query = process_query(query, args.model_path, args.batch_key, seed=args.seed)

    query.write_h5ad(f"{args.study_key}_processed.h5ad")
    print(
        f"Wrote {query.n_obs} cells → "
        f"{args.study_key}_processed.h5ad  (obsm['scvi'] shape: {query.obsm['scvi'].shape})"
    )


if __name__ == "__main__":
    main()
