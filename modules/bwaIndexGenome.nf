/*
 * Create BWA index from a reference genome
 * For BwaAln and BwaMem alignment modules 
 */
process bwaIndexGenome {

    tag "${genomeFasta}"

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }

    container 'variantvalidator/indexgenome:1.1.0'

    // Published index files to directory specified in config: 
    publishDir(path: "$params.resources_dir/indexes/bwa", mode: 'copy')

    input:
    path genomeFasta

    output:
    tuple path(genomeFasta), path("${genomeFasta}.*")

    script:
    """
    echo "Running BWA Indexing"

    # Generate BWA index
    bwa index "${genomeFasta}"

    echo "BWA Genome Indexing complete."
    """
}
