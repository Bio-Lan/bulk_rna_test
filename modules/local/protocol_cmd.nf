process PROTOCOL_CMD {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::python==3.12'
    container "biocontainers/python:3.12"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path assets_dir
    val protocol

    output:
    tuple val(meta), path("${meta.id}.protocol_cmd.txt"), emit: protocol_cmd
    tuple val(meta), path('*.json'),  emit: json

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    def pattern = params.pattern ? "--pattern ${params.pattern}" : ""
    def whitelist = params.whitelist ? "--whitelist \'${params.whitelist}\'" : ""
    """
    protocol_cmd.py \\
        --sample ${prefix} \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} \\
        --well ${params.well} \\
        $pattern \\
        $whitelist
    """
}