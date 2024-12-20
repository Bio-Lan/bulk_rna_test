process STARSOLO_SUMMARY {
    tag "$meta.id"
    label 'process_medium'

    conda 'conda-forge::pandas==1.5.2 bioconda::pysam==0.22.1'
    container "qaqlans/sgr-accura-3"

    input:
    tuple val(meta), path(read_stats), path(summary), path(bam)

    output:
    tuple val(meta), path("multiqc_data/*.json"), emit: json
    tuple val(meta), path("*.counts.txt"), emit: raw_count
    tuple val(meta), path("*.counts_report.txt"), emit: filter_count

    script:
    """
    starsolo_summary.py \\
        --read_stats ${read_stats} \\
        --summary ${summary} \\
        --sample ${meta.id} \\
        --bam ${bam} \\
        --umi_cutoff ${params.umi_cutoff} \\
        --read_cutoff ${params.read_cutoff} \\
        --gene_cutoff ${params.gene_cutoff}
    """
}