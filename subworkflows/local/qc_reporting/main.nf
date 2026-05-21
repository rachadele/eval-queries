/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC_REPORTING subworkflow — query-query edition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - Per-pair MultiQC reports with actual vs. predicted composition bar plots
    - One summary heatmap figure per taxonomy level (F1/precision/recall × scvi/seurat)
----------------------------------------------------------------------------------------
*/

include { PLOT_PREDICTED_COMPOSITION } from "$projectDir/modules/local/plot_predicted_composition/main"
include { PLOT_F1_HEATMAP            } from "$projectDir/modules/local/plot_f1_heatmap/main"
include { MULTIQC                    } from "$projectDir/modules/local/multiqc/main"

workflow QC_REPORTING {
    take:
    predicted_meta_channel  // (method, query_relabel_tsv, ref_path, predicted_meta_tsv)
    f1_scores_channel       // (method, summary.scores.tsv)

    main:
    // ── Per-pair composition bar plots ──────────────────────────────────────
    plot_input = predicted_meta_channel.map { method, query_path, ref_path, predicted_meta ->
        def query_name = query_path.getName().replace('.obs.relabel.tsv', '')
        def ref_name   = ref_path.getName()
                            .replace('_processed.h5ad', '')
                            .replace('_processed.rds', '')
        def study_name = query_name.tokenize('_')[0]
        [method, study_name, query_name, ref_name, predicted_meta]
    }

    PLOT_PREDICTED_COMPOSITION(plot_input)

    // One MultiQC report per (method, study_name, ref_name)
    multiqc_input = PLOT_PREDICTED_COMPOSITION.out.mqc_dir
        .groupTuple(by: [0, 1, 2])

    MULTIQC(multiqc_input)

    // ── Summary heatmaps across all pairs ───────────────────────────────────
    scvi_scores = f1_scores_channel
        .filter { method, scores -> method == "scvi" }
        .map    { method, scores -> scores }
        .collect()

    seurat_scores = f1_scores_channel
        .filter { method, scores -> method == "seurat" }
        .map    { method, scores -> scores }
        .collect()

    PLOT_F1_HEATMAP(scvi_scores, seurat_scores)

    emit:
    multiqc_html  = MULTIQC.out.multiqc_html
    heatmap_pngs  = PLOT_F1_HEATMAP.out.heatmap_pngs
}
