/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEURAT_PIPELINE subworkflow — query-query edition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Takes pre-computed pairwise (query_rds, ref_rds) combos and runs Seurat label
    transfer + classification.
----------------------------------------------------------------------------------------
*/

include { PREDICT_SEURAT } from "$projectDir/modules/local/predict_seurat/main"
include { CLASSIFY_ALL   } from "$projectDir/modules/local/classify_all/main"

workflow SEURAT_PIPELINE {
    take:
    pairwise_combos_seurat  // channel of tuples: (query_rds, ref_rds)
    ref_keys
    ref_region_mapping

    main:
    // Run Seurat label transfer
    PREDICT_SEURAT(pairwise_combos_seurat, ref_keys)

    // Prepare channel for classification
    seurat_scores_channel = PREDICT_SEURAT.out.pred_scores_channel.map { query_path, ref_path, scores_path ->
        ["seurat", query_path, ref_path, scores_path]
    }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, seurat_scores_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
