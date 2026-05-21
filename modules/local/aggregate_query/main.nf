process AGGREGATE_QUERY {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    tuple val(study_key), val(relabel_path), path(h5ad_files), val(batch_key)
    val ref_keys

    output:
    tuple val(study_key), val(relabel_path), path("${study_key}.h5ad"), val(batch_key), emit: aggregated_adata

    script:
    """
    python $projectDir/bin/aggregate_query.py \\
        --study_key ${study_key} \\
        --h5ad_files ${h5ad_files.join(' ')} \\
        --relabel_path ${relabel_path} \\
        --gene_mapping ${params.gene_mapping} \\
        --mapping_file ${params.relabel_r} \\
        --ref_keys ${ref_keys} \\
        --subsample_ref ${params.subsample_ref} \\
        --seed ${params.seed} \\
        ${params.remove_unknown ? '--remove_unknown' : ''}
    """
}
