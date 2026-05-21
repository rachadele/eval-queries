#!/usr/bin/env python3

"""
For one (query, ref) label-transfer result, generate MultiQC-ready TSV files
showing actual vs. predicted cell type counts at each level of granularity.

One TSV is written per ref_key level, named:
    <query_name>_<ref_name>_<level>_actual_vs_predicted_mqc.tsv

Rows = true cell type, columns = predicted cell type.
"""

import argparse
import os
import sys

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate MultiQC bargraph TSVs from a predicted_meta file."
    )
    parser.add_argument("--predicted_meta", type=str, required=True,
                        help="Path to the predictions TSV from classify_all.py")
    parser.add_argument("--query_name", type=str, required=True)
    parser.add_argument("--ref_name", type=str, required=True)
    parser.add_argument("--method", type=str, required=True)
    parser.add_argument("--ref_keys", type=str, nargs="+",
                        default=["subclass", "class", "family", "global"])
    parser.add_argument("--outdir", type=str, default=".")
    return parser.parse_args()


def main():
    args = parse_arguments()

    df = pd.read_csv(args.predicted_meta, sep="\t")

    os.makedirs(args.outdir, exist_ok=True)

    for level in args.ref_keys:
        pred_col = f"predicted_{level}"
        if level not in df.columns or pred_col not in df.columns:
            print(f"Skipping level '{level}': columns not found.", file=sys.stderr)
            continue

        counts = (
            df.groupby([level, pred_col])
            .size()
            .unstack(fill_value=0)
        )
        counts.index.name = level

        out_path = os.path.join(
            args.outdir,
            f"{args.method}_{args.query_name}_{args.ref_name}_{level}_actual_vs_predicted_mqc.tsv",
        )
        counts.to_csv(out_path, sep="\t")
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
