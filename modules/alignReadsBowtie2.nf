/*
 * Align reads to the indexed genome using Bowtie2
 * Bowtie2 is flexible for both short and long reads
 */

 process alignReadsBowtie2 {
    
    tag "${sample_id}"

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'community.wave.seqera.io/library/bowtie2_samtools:3daa9c846bba9a07'

    publishDir(path: "${params.outdir}/alignment/${params.aligner}", mode: 'copy')

    input: 
    tuple val(sample_id), path(reads) // paired-end reads
    path requiredIndexFiles // all .bt2 files from bowtie2-build

    output: 
    tuple val(sample_id), file("${sample_id}.bam")

    script: 

    def read1 = reads[0]
    def read2 = reads.size() > 1 ? reads[1] : null 

    // Pick the first .1.bt2 file from the index directory
    def indexPrefix = requiredIndexFiles.find { it.name.endsWith('.1.bt2') && !it.name.contains('.rev') }.name.replaceAll(/\.1\.bt2$/, '')

    """
    echo "Running Bowtie2 for sample: ${sample_id}"
    echo "Using index prefix: ${indexPrefix}" 

    if [ -f "${read1}" ]; then
        if [ -f "${read2}" ]; then
            # Paired-end mode
            bowtie2 -x ${indexPrefix} -1 ${read1} -2 ${read2} -p ${task.cpus} |
            samtools view -b - |
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
        else
            # Single-end mode
            bowtie2 -x ${indexPrefix} -U ${read1} -p ${task.cpus} |
            samtools view -b - |
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
        fi
    else
        echo "Error: Read file ${read1} does not exist for ${sample_id}."
        exit 1
    fi

    echo "Alignment complete"
    """
 }