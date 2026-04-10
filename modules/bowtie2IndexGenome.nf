/*
 * Create Bowtie2 genome index from FASTA reference
 */

 process bowtie2IndexGenome {

    tag "${genome_file}"

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }

    container 'community.wave.seqera.io/library/bowtie2_samtools:3daa9c846bba9a07'

    input:
    path genome_file

    output:
    path("*.bt2")

    // Published index files to directory specified in config: 
    publishDir(path: "$params.resources_dir/indexes/bowtie2", mode: 'copy')

    script:
    """
    echo "Running Bowtie2 index on genome: ${genome_file}"

    # Extract clean name for genome file 
    GENOME_BASE=\$(basename ${genome_file} .fasta)

    echo "Using basename: \$GENOME_BASE"

    # Create Bowtie2 index
    bowtie2-build ${genome_file} \$GENOME_BASE

    echo "Bowtie2 indexing complete." 
    """
 }