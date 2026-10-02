# Pairwise Query-Query Label Transfer Pipeline

Benchmarks cross-study cell type annotation concordance by using each mouse scRNA-seq study as a reference to annotate all other studies. This is a companion to the [`nextflow_eval_pipeline`](../annotation-benchmark/) (Census reference → query), designed to determine whether annotation failures (e.g. for microglia, OPCs) originate from the reference or from the query data itself.

**Hypothesis:** if query→query label transfer performs well, the problem lies with the Census reference. If it also fails pairwise, the issue is with query data concordance with established cell type signatures.

---

## How it works

1. **Aggregate** all per-sample h5ads for each study into a single study-level h5ad (gene-mapped, relabeled, subsampled to ≤500 cells per cell type).
2. **Project** each study-level h5ad through the Census SCVI model into a shared latent space.
3. **Generate all pairwise combinations** (N studies → N×(N−1) directional pairs, excluding self-pairs).
4. **Transfer labels** using two independent methods in parallel:
   - **SCVI + Random Forest**: trains an RF classifier on the reference SCVI embeddings, predicts cell types in the query.
   - **Seurat**: performs anchor-based label transfer (FindTransferAnchors + TransferData).
5. **Classify and evaluate** predictions at four levels of granularity (`subclass`, `class`, `family`, `global`) using F1, accuracy, NMI, ARI, and confusion matrices.

---

## Studies

The 6 relabeled mouse studies (30 directional pairs):

| Study | Brain Region | Relabel file |
|-------|-------------|--------------|
| GSE124952 | prefrontal cortex | `GSE124952_relabel.tsv` |
| GSE181021.2 | prefrontal cortex | `GSE181021.2_relabel.tsv` |
| GSE185454 | hippocampus | `GSE185454_relabel.tsv` |
| GSE214244.1 | entorhinal cortex | `GSE214244.1_relabel.tsv` |
| GSE247339.1 | hippocampus | `GSE247339.1_relabel.tsv` |
| GSE247339.2 | cortex | `GSE247339.2_relabel.tsv` |

---

## Pipeline structure

```
eval-queries/
├── main.nf                                  # Main workflow
├── nextflow.config                          # Pipeline configuration
├── conf/
│   ├── base.config -> (symlink)             # SLURM executor settings
│   ├── params.config                        # Default parameters
│   └── modules.config                       # Module-level config
├── modules/local/
│   ├── aggregate_query/main.nf              # Concat + relabel + subsample per study
│   ├── map_aggregated_query/main.nf         # QC + SCVI projection
│   ├── classify_all/main.nf                 # Metrics + confusion matrices
│   ├── run_setup/ -> (symlink)              # Download Census SCVI model
│   ├── scvi_predict/ -> (symlink)             # Random Forest prediction
│   ├── predict_seurat/ -> (symlink)         # Seurat label transfer
│   └── query_process_seurat/ -> (symlink)   # Convert h5ad → Seurat RDS
├── subworkflows/local/
│   ├── scvi_pipeline/main.nf               # SCVI RF pipeline (pairwise combos)
│   └── seurat_pipeline/main.nf             # Seurat pipeline (pairwise combos)
├── bin/
│   ├── aggregate_query.py                   # Study aggregation script
│   ├── map_aggregated_query.py              # QC + SCVI projection script
│   ├── utils.py -> (symlink)               # Shared utilities
│   ├── predict_scvi.py -> (symlink)
│   ├── predict_seurat.R -> (symlink)
│   ├── classify_all.py -> (symlink)
│   ├── seurat_preprocessing.R -> (symlink)
│   └── seurat_functions.R -> (symlink)
├── assets/
│   └── gemma_genes.tsv -> (symlink)        # ENSEMBL → gene symbol mapping
└── meta/
    ├── query_region_mapping.yaml            # Study → brain region
    └── census_map_mouse_author.tsv -> (symlink)  # Label hierarchy
```

---

## Input data

| Parameter | Description | Example path |
|-----------|-------------|--------------|
| `queries_adata` | Glob for per-sample h5ad files | `/path/to/h5ad/**/*.h5ad` |
| `relabel_q` | Glob for `<study_key>_relabel.tsv` files | `/path/to/relabel_mus_musculus/*_relabel.tsv` |
| `relabel_r` | Label hierarchy mapping TSV | `$projectDir/meta/census_map_mouse_author.tsv` |

**File naming convention for h5ads:** `<study_key>_<anything>.h5ad`
(e.g. `GSE247339.2_1051970_GSM7887408.h5ad`). The pipeline extracts the study key from the first `_`-delimited token.

---

## Running the pipeline

### Prerequisites

- Nextflow ≥ 23.04
- Conda environments:
  - `scanpyenv` — Python (scanpy, scvi-tools, sklearn, anndata)
  - `r4.3` — R (Seurat, SeuratDisk)
- SLURM cluster access (configured in `conf/base.config`)

### Create a params file

```json
{
  "queries_adata": "/space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true/h5ad/**/*.h5ad",
  "relabel_q": "/space/grp/rschwartz/rschwartz/annotation-benchmark/meta/relabel_mus_musculus/*_relabel.tsv",
  "relabel_r": "/space/grp/rschwartz/rschwartz/annotation-benchmark/meta/census_map_mouse_author.tsv"
}
```

### Run

```bash
nextflow run main.nf \
  -params-file params.mm.json \
  -profile conda \
  -resume
```

### Smoke test (2 studies, small subsample)

```bash
nextflow run main.nf \
  -params-file params.mm.json \
  -profile conda \
  --subsample_ref 50 \
  -resume
```

---

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `organism` | `mus_musculus` | Species for Census SCVI model |
| `census_version` | `2025-01-30` | Census release to pull the model from |
| `subsample_ref` | `500` | Max cells per cell type when building each study reference |
| `ref_keys` | `[subclass, class, family, global]` | Label granularity levels to evaluate |
| `cutoff` | `0` | Minimum confidence score to assign a label (0 = no threshold) |
| `use_gap` | `true` | Use gap-based confidence classification instead of raw probability cutoff |
| `nmads` | `5` | MADs for QC outlier detection |
| `remove_unknown` | `true` | Drop cells with no label in the relabel file |
| `normalization_method` | `SCT` | Seurat normalisation method |
| `seed` | `42` | Random seed |

Batch keys per study are configured in `conf/params.config` under `batch_keys`.

---

## Output structure

```
query_query/mus_musculus/<subsample_ref>/SCT/gap_<use_gap>/cutoff_<cutoff>/
├── scvi/
│   └── <query_study>/
│       └── <ref_study>/
│           └── <query_name>/
│               ├── label_transfer_metrics/*.summary.scores.tsv
│               ├── confusion/*.png
│               └── predicted_meta/*.predictions.*.tsv
└── seurat/
    └── <query_study>/
        └── <ref_study>/
            └── <query_name>/
                └── [same structure]
```

With 6 studies, expect **30 result directories per method** (60 total).

---

## Relationship to `nextflow_eval_pipeline`

This pipeline reuses the same SCVI model, bin scripts, and Nextflow modules as the original Census reference → query pipeline. The key differences are:

| | `nextflow_eval_pipeline` | `eval-queries` |
|---|---|---|
| Reference | CellXGene Census datasets | Query studies (aggregated) |
| Query | Per-sample h5ads | Per-study aggregated h5ads |
| Pairs | Census ref × each query | All pairwise query × query |
| Aggregation step | None (Census handled upstream) | `aggregate_query.py` |
| SCVI projection | `process_query.py` (includes relabeling) | `map_aggregated_query.py` (projection only) |
