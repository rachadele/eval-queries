process MULTIQC {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/${method}/${study_name}/${ref_name}", mode: 'copy', pattern: "*multiqc_report.html"

    input:
    tuple val(method), val(study_name), val(ref_name), path(mqc_dirs)

    output:
    tuple val(method), val(study_name), val(ref_name), path("*multiqc_report.html"), emit: multiqc_html

    script:
    """
    cp ${params.multiqc_config} new_config.yaml
    echo 'title: "${study_name} vs ${ref_name} [${method}]"' >> new_config.yaml
    multiqc . --config new_config.yaml
    """
}
