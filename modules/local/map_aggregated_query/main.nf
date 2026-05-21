process MAP_AGGREGATED_QUERY {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val model_path
    tuple val(study_key), val(relabel_path), path(h5ad_file), val(batch_key)

    output:
    path "${study_key}_processed.h5ad", emit: processed_adata

    script:
    """
    python $projectDir/bin/map_aggregated_query.py \\
        --h5ad_path ${h5ad_file} \\
        --study_key ${study_key} \\
        --model_path ${model_path} \\
        --batch_key ${batch_key} \\
        --seed ${params.seed}
    """
}
