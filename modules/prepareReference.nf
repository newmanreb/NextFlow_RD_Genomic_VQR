/*
 * Prepare reference genome assets for variant calling
 * Generates: .fai and .dict
 */

 process prepareReference {

    tag "${fasta}"

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }

    container 'variantvalidator/indexgenome:1.1.0'

    publishDir "$params.resources_dir/genomes", mode: 'copy'

    input:
    path fasta

    output:
    tuple path(fasta), path("*.fai"), path("*.dict")

    script:
    """
    echo "Preparing reference genome"

    samtools faidx ${fasta}

    BASENAME=\$(basename ${fasta} .fasta)

    picard CreateSequenceDictionary \
        R=${fasta} \
        O=\${BASENAME}.dict

    echo "Reference preparation complete"
    """
}