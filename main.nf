#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pairwise Query-Query Label Transfer Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uses each mouse query study as a reference to annotate all other mouse query
    studies, reusing the same Census SCVI model for projection.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_SETUP            } from "$projectDir/modules/local/run_setup/main"
include { AGGREGATE_QUERY      } from "$projectDir/modules/local/aggregate_query/main"
include { MAP_AGGREGATED_QUERY } from "$projectDir/modules/local/map_aggregated_query/main"
include { QUERY_PROCESS_SEURAT } from "$projectDir/modules/local/query_process_seurat/main"
include { SCVI_PIPELINE        } from "$projectDir/subworkflows/local/scvi_pipeline/main"
include { SEURAT_PIPELINE      } from "$projectDir/subworkflows/local/seurat_pipeline/main"
include { QC_REPORTING         } from "$projectDir/subworkflows/local/qc_reporting/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // 1. Download Census SCVI model
    RUN_SETUP(params.organism, params.census_version)

    // 2. Build (study_key → [h5ad files]) channel
    //    File naming convention: <study_key>_<sample_id>.h5ad
    query_paths = Channel.fromPath(params.queries_adata)
        .map { p ->
            def study_key = p.getName().tokenize('_')[0]
            [study_key, p]
        }
        .groupTuple()   // [ study_key, [h5ad, h5ad, ...] ]

    // Map relabel TSVs: <study_key>_relabel.tsv
    relabel_paths = Channel.fromPath(params.relabel_q)
        .map { p ->
            def study_key = p.getName().replace('_relabel.tsv', '')
            [study_key, p]
        }

    // Join h5ad list with relabel path by study_key, add batch_key
    combined = query_paths
        .combine(relabel_paths, by: 0)
        .map { study_key, h5ad_files, relabel_path ->
            def batch_key = params.batch_keys.get(study_key, "sample_id")
            [study_key, relabel_path, h5ad_files, batch_key]
        }

    // 3. Aggregate all samples per study into one h5ad
    AGGREGATE_QUERY(combined, params.ref_keys.join(' '))

    // 4. QC + SCVI projection
    MAP_AGGREGATED_QUERY(RUN_SETUP.out, AGGREGATE_QUERY.out.aggregated_adata)

    processed_adatas = MAP_AGGREGATED_QUERY.out.processed_adata

    // 5. Generate all pairwise combos (exclude self-pairs) for SCVI
    pairwise_combos_adata = processed_adatas
        .combine(processed_adatas)
        .filter { query_h5ad, ref_h5ad ->
            query_h5ad.getName() != ref_h5ad.getName()
        }

    // 6. SCVI-based RF prediction pipeline
    SCVI_PIPELINE(
        pairwise_combos_adata,
        params.ref_keys.join(' '),
        params.ref_region_mapping
    )

    // 7. Seurat pipeline: convert all processed h5ads to .rds, then pairwise combos
    all_seurat = QUERY_PROCESS_SEURAT(processed_adatas)

    pairwise_combos_seurat = all_seurat
        .combine(all_seurat)
        .filter { query_rds, ref_rds ->
            query_rds.getName() != ref_rds.getName()
        }

    SEURAT_PIPELINE(
        pairwise_combos_seurat,
        params.ref_keys.join(' '),
        params.ref_region_mapping
    )

    // 8. QC reports and summary heatmaps
    predicted_meta_combined = SCVI_PIPELINE.out.predicted_meta.concat(
        SEURAT_PIPELINE.out.predicted_meta
    )

    f1_scores_combined = SCVI_PIPELINE.out.f1_scores.concat(
        SEURAT_PIPELINE.out.f1_scores
    )

    QC_REPORTING(predicted_meta_combined, f1_scores_combined)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    println(workflow.success ?
        """
        ================================================================================
        Pipeline execution summary
        --------------------------------------------------------------------------------

        Run as      : ${workflow.commandLine}
        Started at  : ${workflow.start}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        Config files: ${workflow.configFiles}
        Exit status : ${workflow.exitStatus}
        Output dir  : ${params.outdir}

        --------------------------------------------------------------------------------
        ================================================================================
        """.stripIndent() :
        """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """.stripIndent()
    )
}

workflow.onError = {
    println "Error: something went wrong, check the pipeline log at '.nextflow.log'"
}
