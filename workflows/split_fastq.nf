/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process split_fastq {
    tag "$meta.id"
    label 'process_medium'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 conda-forge::xopen==2.0.1'
    container "qaqlans/sgrdocker_accura_tools1"
    
    input:
    tuple val(meta), path(reads)
    path assets_dir
    val protocol
    path split_inf

    output:
    tuple val(meta), path("${meta.id}/")
    tuple val(meta), path("${meta.id}.*.json")

    script:
    def prefix = "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def args = task.ext.args ?: ''
    """
    split_fastq.py \\
        --sample $prefix \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
        --split_inf $split_inf \\
        --assets_dir $assets_dir \\
        --protocol $protocol \\
        --well ${params.well} \\
        --pattern ${params.pattern} \\
        --whitelist \"${params.whitelist}\" \\
        $args
    """
}

workflow SPLITFASTQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    split_fastq (
        ch_samplesheet,
        "${projectDir}/assets/",
        params.protocol,
        params.split_inf,
    )
}