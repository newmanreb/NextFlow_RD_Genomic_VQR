/*
 * Align reads to the indexed genome. BwaMem is designed for longer reads 
 * (>70-100bp) and is a more modern, standard default. 
 */
process alignReadsBwaMem {

    tag "${sample_id}"

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'community.wave.seqera.io/library/bwa_samtools:3704450416e4d5eb'

    publishDir(path: "${params.outdir}/alignment/bwa_mem", mode: 'copy')

    input:
    tuple val(sample_id), path(reads)   // reads is a tuple of paths for paired-end reads
    path requiredIndexFiles

    output:
    tuple val(sample_id), file("${sample_id}.bam")

    script:
    """
    echo "Running BWA-MEM for sample: ${sample_id}" 

    # Derive index prefix from .amb file 
    INDEX=\$(ls *.amb | head -n 1 | sed 's/\\.amb\$//')
    #INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')

    echo "Using index prefix: \$INDEX"

    # Check if the input FASTQ files exist
    if [ -f "${reads[0]}" ]; then
        if [ -f "${reads[1]}" ]; then
            # Paired-end mode
            bwa mem -M -k 16 -t ${task.cpus} \$INDEX ${reads[0]} ${reads[1]} |
            samtools view -b - |
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
        else
            # Single-end mode
            bwa mem -M -k 16 -t ${task.cpus} \$INDEX ${reads[0]} |
            samtools view -b - |
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
        fi
    else
        echo "Error: Read file ${reads[0]} does not exist for sample ${sample_id}."
        exit 1
    fi

    echo "Alignment complete"
    """
}
