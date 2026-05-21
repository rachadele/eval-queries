process PLOT_PREDICTED_COMPOSITION {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    tuple val(method), val(study_name), val(query_name), val(ref_name), path(predicted_meta)

    output:
    tuple val(method), val(study_name), val(ref_name), path("${method}_${query_name}_${ref_name}/"), emit: mqc_dir

    script:
    """
    python $projectDir/bin/plot_predicted_composition.py \\
        --predicted_meta ${predicted_meta} \\
        --query_name ${query_name} \\
        --ref_name ${ref_name} \\
        --method ${method} \\
        --ref_keys ${params.ref_keys.join(' ')} \\
        --outdir ${method}_${query_name}_${ref_name}
    """
}
