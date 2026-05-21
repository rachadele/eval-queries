process PLOT_F1_HEATMAP {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path scvi_scores,   stageAs: 'scvi/*'
    path seurat_scores, stageAs: 'seurat/*'

    output:
    path "*.png", emit: heatmap_pngs
    path "*.tsv", emit: heatmap_tsvs

    script:
    """
    python $projectDir/bin/plot_f1_heatmap.py \\
        --scvi_scores ${scvi_scores} \\
        --seurat_scores ${seurat_scores} \\
        --ref_keys ${params.ref_keys.collect { "\"${it}\"" }.join(' ')} \\
        --target_celltypes ${params.target_celltypes.collect { "\"${it}\"" }.join(' ')}
    """
}
