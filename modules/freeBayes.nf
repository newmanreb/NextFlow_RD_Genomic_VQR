process freeBayes {
    
    container 'community.wave.seqera.io/library/freebayes:1.3.10--a4e01a89d7090203'

    tag "$sample_id"

    publishDir(path: "${params.outdir}/VCF", mode: 'copy')

    input: 
    tuple val(sample_id), file(bam), file(bai)
    file fasta

    output: 
    tuple val(sample_id), file("${sample_id}.vcf")

    script: 
    """
    freebayes \
        -f ${fasta} \
        ${bam} \
        > ${sample_id}.vcf
    """
}