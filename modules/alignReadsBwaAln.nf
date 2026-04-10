/*
 * Align reads to the indexed genome. BwaAln is designed for shorter reads 
 * (35-75 bp) and is an older, largely legacy tool. 
 */

process alignReadsBwaAln {

    tag "${sample_id}"

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'community.wave.seqera.io/library/bwa_samtools:3704450416e4d5eb'

    publishDir(path: "${params.outdir}/alignment/${params.aligner}", mode: 'copy')

    input:
    tuple val(sample_id), path(reads)   // reads is a tuple of paths for paired-end reads
    path requiredIndexFiles

    output:
    tuple val(sample_id), file("${sample_id}.bam")

    script:

    def read1 = reads[0]
    def read2 = reads.size() > 1 ? reads[1] : null

    """
    echo "Running BWA-ALN for sample: ${sample_id}"

    # Find the BWA index
    INDEX=\$(ls *.amb | head -n 1 | sed 's/\\.amb\$//')
    echo "Using index prefix: \$INDEX"

    if [ -f "${read1}" ]; then
        if [ -f "${read2}" ]; then 
            echo "Paired-end mode" 

            bwa aln -t ${task.cpus} \$INDEX ${read1} > ${sample_id}_1.sai 
            bwa aln -t ${task.cpus} \$INDEX ${read2} > ${sample_id}_2.sai 

            bwa sampe \$INDEX \
                ${sample_id}_1.sai ${sample_id}_2.sai \
                ${read1} ${read2} |
            samtools view -b - |
            samtools addreplacerg \
                -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" \
                - > ${sample_id}.bam
        else 
            echo "Single-end mode" 

            bwa aln -t ${task.cpus} \$INDEX ${read1} > ${sample_id}.sai 

            bwa samse \$INDEX ${sample_id}.sai ${read1} |
            samtools view -b - |
            samtools addreplacerg \
                -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" \
                - > ${sample_id}.bam
        fi 
    
    else 
        echo "Error: Read file ${read1} does not exist for ${sample_id}."
        exit 1 
    fi

    echo "Alignment complete"
    """
}
