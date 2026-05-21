# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

Pairwise query-query label transfer benchmarking pipeline. Each of 7 relabeled mouse scRNA-seq studies is used as a reference to annotate all other studies (42 directional pairs). This is a companion to `../nextflow_eval_pipeline/` (Census → query), designed to isolate whether annotation failures originate from the reference or the query data.

## Running the pipeline

```bash
# Standard run
nextflow run main.nf -params-file params.mm.json -profile conda -resume

# Smoke test (fast, small subsample)
nextflow run main.nf -params-file params.mm.json -profile conda --subsample_ref 50 -resume
```

Requires SLURM cluster access. Work directory: `/cosmos/data/nextflow-eval-pipeline-work`.

## Architecture

### Data flow

```
per-sample h5ads (grouped by study_key)
  → AGGREGATE_QUERY        # concat + relabel + subsample ≤500 cells/type
  → MAP_AGGREGATED_QUERY   # QC + SCVI projection (obsm['scvi'])
  → pairwise combos (N×(N−1), no self-pairs)
  → SCVI_PIPELINE          # RF classifier on SCVI embeddings
  → SEURAT_PIPELINE        # Seurat anchor-based label transfer
  → QC_REPORTING           # composition plots + heatmaps
```

### File naming convention

Input h5ads must be named `<study_key>_<anything>.h5ad`. The pipeline extracts `study_key` from the first `_`-delimited token (e.g., `GSE247339.2_1051970_GSM7887408.h5ad` → `GSE247339.2`).

Aggregated files: `<study_key>.h5ad` → `<study_key>_processed.h5ad` after SCVI. `classify_all/main.nf` strips `_processed.h5ad` to recover the ref_name.

### Label hierarchy

Four levels of granularity evaluated: `subclass → class → family → global`. Mapping defined in `meta/census_map_mouse_author.tsv`. Each study has a study-specific relabel TSV under `meta/relabel_mus_musculus/`.

### Symlinked components

Many components are symlinked from `../nextflow_eval_pipeline/`:
- `bin/`: `utils.py`, `predict_scvi.py`, `predict_seurat.R`, `classify_all.py`, `setup.py`, seurat R scripts
- `modules/local/`: `run_setup`, `rf_predict`, `predict_seurat`, `query_process_seurat`
- `assets/gemma_genes.tsv`, `meta/census_map_mouse_author.tsv`, `conf/base.config`

Modifying these symlinked files will affect the original pipeline too. Only `aggregate_query.py`, `map_aggregated_query.py`, `plot_f1_heatmap.py`, `plot_predicted_composition.py`, and their corresponding modules/subworkflows are unique to this pipeline.

## Key parameters (conf/params.config)

| Parameter | Default | Notes |
|-----------|---------|-------|
| `subsample_ref` | `500` | Max cells per cell type per study |
| `ref_keys` | `[subclass, class, family, global]` | Levels to evaluate |
| `use_gap` | `true` | Gap-based confidence vs. raw probability cutoff |
| `cutoff` | `0` | Min confidence threshold |
| `batch_keys` | per-study map | Study-specific AnnData obs column for batch |
| `census_version` | `2025-01-30` | Census SCVI model version |

Batch keys per study are hardcoded in `conf/params.config` and must be updated when adding new studies.

## Conda environments

- Python: `/home/rschwartz/anaconda3/envs/scanpyenv`
- R: `/home/rschwartz/anaconda3/envs/r4.3`

## Output structure

```
query_query/mus_musculus/<subsample_ref>/SCT/gap_<use_gap>/cutoff_<cutoff>/
├── scvi/<query>/<ref>/<query>/
│   ├── label_transfer_metrics/*.summary.scores.tsv
│   ├── confusion/*.png
│   └── predicted_meta/*.predictions.*.tsv
├── seurat/  (same structure)
└── figures/
    ├── {subclass,class,family,global}_metrics.png  # per-level summary heatmaps
    └── pairwise_{f1_score,precision,recall}.png    # target cell type detail
```

## Adding a new study

1. Add study-specific relabel TSV to `meta/relabel_mus_musculus/<study_key>_relabel.tsv`
2. Add `<study_key>` → brain region entry in `meta/query_region_mapping.yaml`
3. Add `<study_key>` → batch_key entry in `conf/params.config` under `params.batch_keys`
4. Ensure h5ad files follow the `<study_key>_<sample_id>.h5ad` naming convention
5. Update `params.mm.json` if the glob needs to capture new study files
