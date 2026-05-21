#!/usr/bin/env python3

"""
Generate two sets of heatmap figures for query-query label transfer results:

1. Per-level summary (averaged across reference studies)
   One PNG per taxonomy level, 3 × 2 grid (f1|precision|recall × scvi|seurat).
   Rows = cell type labels; columns = query study.

2. Pairwise detail for target cell types (Microglia, OPC by default)
   One PNG per metric, 1 × 2 grid (scvi|seurat).
   Rows = query → ref pairs; columns = target cell type.
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


METRICS = ["f1_score", "precision", "recall"]
METRIC_LABELS = {"f1_score": "F1", "precision": "Precision", "recall": "Recall"}
METHODS = ["scvi", "seurat"]


def load_scores(files, method_tag):
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            df["method"] = method_tag
            dfs.append(df)
        except Exception as e:
            print(f"Warning: could not read {f}: {e}", file=sys.stderr)
    if not dfs:
        return pd.DataFrame()
    combined = pd.concat(dfs, ignore_index=True)
    zero_support = combined["support"] == 0
    combined.loc[zero_support, ["f1_score", "recall"]] = float("nan")
    return combined


# ── Plot 1: per-level summary averaged across refs ────────────────────────────

def make_summary_pivot(df, level, method, metric):
    sub = df[(df["key"] == level) & (df["method"] == method) & (df["support"] > 0)].copy()
    if sub.empty:
        return pd.DataFrame()
    return (
        sub.groupby(["label", "study"])[metric]
        .mean()
        .unstack("study")
        .sort_index()
    )


def save_summary_tsv(df, level):
    rows = []
    for method in METHODS:
        sub = df[(df["key"] == level) & (df["method"] == method) & (df["support"] > 0)].copy()
        if sub.empty:
            continue
        for metric in METRICS:
            stats = (
                sub.groupby(["label", "study"])[metric]
                .agg(mean="mean", sd="std")
                .reset_index()
            )
            stats["method"] = method
            stats["metric"] = metric
            rows.append(stats)
    if rows:
        result = pd.concat(rows, ignore_index=True)
        result = result[["method", "metric", "label", "study", "mean", "sd"]]
        out_path = f"{level}_metrics.tsv"
        result.to_csv(out_path, sep="\t", index=False)
        print(f"Wrote {out_path}")


def plot_summary(df, ref_keys):
    for level in ref_keys:
        fig, axes = plt.subplots(
            nrows=len(METRICS),
            ncols=len(METHODS),
            figsize=(7 * len(METHODS), 4 * len(METRICS)),
        )
        fig.suptitle(f"Label transfer metrics — {level}", fontsize=14, y=1.01)

        for col, method in enumerate(METHODS):
            for row, metric in enumerate(METRICS):
                ax = axes[row][col]
                pivot = make_summary_pivot(df, level, method, metric)

                if pivot.empty:
                    ax.set_visible(False)
                    continue

                sns.heatmap(
                    pivot,
                    ax=ax,
                    vmin=0, vmax=1,
                    cmap="RdYlGn",
                    annot=True,
                    fmt=".2f",
                    linewidths=0.4,
                    cbar_kws={"shrink": 0.7},
                )
                ax.set_title(f"{method.upper()} — {METRIC_LABELS[metric]}", fontsize=11)
                ax.set_xlabel("Query study")
                ax.set_ylabel("Cell type")
                ax.tick_params(axis="x", rotation=45)
                ax.tick_params(axis="y", rotation=0)

        plt.tight_layout()
        out_path = f"{level}_metrics.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_path}")
        save_summary_tsv(df, level)


# ── Plot 2: per-level reference summary (ref × query, averaged across cell types) ─

def make_ref_pivot(df, level, method, metric):
    sub = df[(df["key"] == level) & (df["method"] == method) & (df["support"] > 0)].copy()
    if sub.empty:
        return pd.DataFrame()
    return (
        sub.groupby(["reference", "study"])[metric]
        .mean()
        .unstack("study")
        .sort_index()
    )


def save_ref_summary_tsv(df, level):
    rows = []
    for method in METHODS:
        sub = df[(df["key"] == level) & (df["method"] == method) & (df["support"] > 0)].copy()
        if sub.empty:
            continue
        for metric in METRICS:
            stats = (
                sub.groupby(["reference", "study"])[metric]
                .agg(mean="mean", sd="std")
                .reset_index()
            )
            stats["method"] = method
            stats["metric"] = metric
            rows.append(stats)
    if rows:
        result = pd.concat(rows, ignore_index=True)
        result = result[["method", "metric", "reference", "study", "mean", "sd"]]
        out_path = f"{level}_ref_summary.tsv"
        result.to_csv(out_path, sep="\t", index=False)
        print(f"Wrote {out_path}")


def plot_ref_summary(df, ref_keys):
    for level in ref_keys:
        fig, axes = plt.subplots(
            nrows=len(METRICS),
            ncols=len(METHODS),
            figsize=(7 * len(METHODS), 4 * len(METRICS)),
        )
        fig.suptitle(f"Reference performance — {level}", fontsize=14, y=1.01)

        for col, method in enumerate(METHODS):
            for row, metric in enumerate(METRICS):
                ax = axes[row][col]
                pivot = make_ref_pivot(df, level, method, metric)

                if pivot.empty:
                    ax.set_visible(False)
                    continue

                sns.heatmap(
                    pivot,
                    ax=ax,
                    vmin=0, vmax=1,
                    cmap="RdYlGn",
                    annot=True,
                    fmt=".2f",
                    linewidths=0.4,
                    cbar_kws={"shrink": 0.7},
                )
                ax.set_title(f"{method.upper()} — {METRIC_LABELS[metric]}", fontsize=11)
                ax.set_xlabel("Query study")
                ax.set_ylabel("Reference study")
                ax.tick_params(axis="x", rotation=45)
                ax.tick_params(axis="y", rotation=0)

        plt.tight_layout()
        out_path = f"{level}_ref_summary.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_path}")
        save_ref_summary_tsv(df, level)


# ── Plot 3: pairwise detail for target cell types ─────────────────────────────

def make_pairwise_pivot(df, method, metric, target_labels):
    method_df = df[df["method"] == method].copy()
    method_df["pair"] = method_df["reference"] + " → " + method_df["study"]
    all_pairs = sorted(method_df["pair"].unique())

    sub = method_df[
        (method_df["key"] == "subclass") &
        (method_df["label"].isin(target_labels)) &
        (method_df["support"] > 0)
    ]

    if sub.empty:
        return pd.DataFrame(index=all_pairs, columns=target_labels)

    return (
        sub.pivot_table(index="pair", columns="label", values=metric, aggfunc="mean")
        .reindex(index=all_pairs, columns=target_labels)
        .sort_index()
    )


def save_pairwise_tsv(df, metric, target_labels):
    rows = []
    for method in METHODS:
        method_df = df[df["method"] == method].copy()
        method_df["pair"] = method_df["reference"] + " → " + method_df["study"]
        sub = method_df[
            (method_df["key"] == "subclass") &
            (method_df["label"].isin(target_labels)) &
            (method_df["support"] > 0)
        ]
        if sub.empty:
            continue
        stats = (
            sub.groupby(["pair", "label"])[metric]
            .agg(mean="mean", sd="std")
            .reset_index()
        )
        stats["method"] = method
        rows.append(stats)
    if rows:
        result = pd.concat(rows, ignore_index=True)
        result = result[["method", "pair", "label", "mean", "sd"]]
        out_path = f"pairwise_{metric}.tsv"
        result.to_csv(out_path, sep="\t", index=False)
        print(f"Wrote {out_path}")


def plot_pairwise(df, target_labels):
    for metric in METRICS:
        fig, axes = plt.subplots(
            nrows=1,
            ncols=len(METHODS),
            figsize=(5 * len(METHODS), 0.4 * 30 + 2),  # scale rows to n pairs
        )
        fig.suptitle(
            f"Pairwise label transfer — {METRIC_LABELS[metric]}: "
            f"{', '.join(target_labels)}",
            fontsize=13,
        )

        for col, method in enumerate(METHODS):
            ax = axes[col]
            pivot = make_pairwise_pivot(df, method, metric, target_labels)

            if pivot.empty:
                ax.set_visible(False)
                continue

            sns.heatmap(
                pivot,
                ax=ax,
                vmin=0, vmax=1,
                cmap="RdYlGn",
                annot=True,
                fmt=".2f",
                linewidths=0.5,
                cbar=False,
            )
            ax.set_title(method.upper(), fontsize=11)
            ax.set_xlabel("Cell type")
            ax.tick_params(axis="x", rotation=30)
            ax.tick_params(axis="y", rotation=0)

            if col == 0:
                ax.set_ylabel("Reference → Query")
            else:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)

        # Shared colorbar on the right
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
        sm = plt.cm.ScalarMappable(cmap="RdYlGn", norm=plt.Normalize(vmin=0, vmax=1))
        fig.colorbar(sm, cax=cbar_ax)
        out_path = f"pairwise_{metric}.png"
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_path}")
        save_pairwise_tsv(df, metric, target_labels)


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate F1/precision/recall heatmaps for query-query label transfer."
    )
    parser.add_argument("--scvi_scores",   type=str, nargs="+", required=True)
    parser.add_argument("--seurat_scores", type=str, nargs="+", required=True)
    parser.add_argument(
        "--ref_keys", type=str, nargs="+",
        default=["subclass", "class", "family", "global"],
    )
    parser.add_argument(
        "--target_celltypes", type=str, nargs="+",
        default=["Microglia", "OPC"],
        help="Cell types for the pairwise detail plot.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    scvi_df   = load_scores(args.scvi_scores,   "scvi")
    seurat_df = load_scores(args.seurat_scores, "seurat")
    df = pd.concat([scvi_df, seurat_df], ignore_index=True)

    if df.empty:
        print("No data found.", file=sys.stderr)
        sys.exit(1)

    plot_summary(df, args.ref_keys)
    plot_ref_summary(df, args.ref_keys)
    plot_pairwise(df, args.target_celltypes)


if __name__ == "__main__":
    main()
