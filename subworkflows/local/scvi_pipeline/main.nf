/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCVI_PIPELINE subworkflow — query-query edition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Takes pre-computed pairwise (query_h5ad, ref_h5ad) combos and runs RF prediction
    + classification.
----------------------------------------------------------------------------------------
*/

include { SCVI_PREDICT } from "$projectDir/modules/local/rf_predict/main"
include { CLASSIFY_ALL } from "$projectDir/modules/local/classify_all/main"

workflow SCVI_PIPELINE {
    take:
    pairwise_combos_adata   // channel of tuples: (query_h5ad, ref_h5ad)
    ref_keys
    ref_region_mapping

    main:
    // Run RF prediction on SCVI embeddings
    SCVI_PREDICT(pairwise_combos_adata, ref_keys)

    // SCVI_PREDICT emits (obs.relabel.tsv, ref_path, rf probs, knn probs).
    // Score each classifier separately, as annotation-benchmark does.
    adata_probs_channel = SCVI_PREDICT.out.probs_channel.flatMap { query_path, ref_path, rf_probs, knn_probs ->
        [
            ["scvi_rf",  query_path, ref_path, rf_probs],
            ["scvi_knn", query_path, ref_path, knn_probs]
        ]
    }

    // Classify and compute metrics
    CLASSIFY_ALL(ref_keys, adata_probs_channel, ref_region_mapping)

    emit:
    f1_scores      = CLASSIFY_ALL.out.f1_score_channel
    predicted_meta = CLASSIFY_ALL.out.predicted_meta_channel
}
