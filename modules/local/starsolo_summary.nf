process STARSOLO_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    input:
    tuple val(meta), path(read_stats), path(summary)
    path assets_dir
    val protocol
    val umi_cutoff
    val read_cutoff
    val gene_cutoff

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.counts.txt"), emit: raw_count
    tuple val(meta), path("*.counts_report.txt"), emit: filter_count

    script:
    """
    starsolo_summary.py \\
        --read_stats ${read_stats} \\
        --summary ${summary} \\
        --sample ${meta.id} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} \\
        --umi_cutoff ${umi_cutoff} \\
        --read_cutoff ${read_cutoff} \\
        --gene_cutoff ${gene_cutoff}
    """
}