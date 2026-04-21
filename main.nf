// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
log.info """\
    ============================================
          DNASeq Pipeline Configuration
    ============================================
    platform        : ${params.platform}
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    qsr truth vcfs  : ${params.qsrVcfs}
    output directory: ${params.outdir}
    fastp           : ${params.fastp}
    fastqc          : ${params.fastqc}
    aligner         : ${params.aligner}
    variant caller  : ${params.variant_caller}
    bqsr            : ${params.bqsr}
    degraded_dna    : ${params.degraded_dna}
    variant_recalibration: ${params.variant_recalibration}
    identity_analysis: ${params.identity_analysis}
    ============================================
""".stripIndent()

// Simple include modules
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
include { prepareReference } from './modules/prepareReference'

// Conditional include modules if selected 
if (params.bqsr) {
    include { baseRecalibrator } from './modules/BQSR'
}
if (params.fastp) {
    include { fastp } from './modules/fastp'
}
if (params.fastqc) {
    include { FASTQC as FASTQC_RAW } from './modules/FASTQC'
    include { FASTQC as FASTQC_TRIMMED } from './modules/FASTQC'
}
if (params.degraded_dna) {
    include { mapDamage2 } from './modules/mapDamage'
    include { indexMapDamageBam } from './modules/indexBam'
}

// Include alignment modules based on chosen aligner 
switch(params.aligner) {
    case 'bwa-mem':
        include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
        include { bwaIndexGenome } from './modules/bwaIndexGenome'
        break
    case 'bwa-aln': 
        include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
        include { bwaIndexGenome } from './modules/bwaIndexGenome'
        break
    case 'bowtie2':
        include { alignReadsBowtie2 } from './modules/alignReadsBowtie2'
        include { bowtie2IndexGenome } from './modules/bowtie2IndexGenome'
        break
    default: 
        error "Unsupported aligner: ${params.aligner}. Please specify 'bwa-mem', 'bwa-aln', or 'bowtie2'."
}

if (params.variant_recalibration) {
    include { variantRecalibrator } from './modules/variantRecalibrator'
} else {
    include { filterVCF } from './modules/filterVCF'
}
if (params.identity_analysis) {
    include { identityAnalysis } from './modules/identityAnalysis'
}

switch(params.variant_caller) {
    case 'haplotype-caller':
        include { haplotypeCaller } from './modules/haplotypeCaller'
        include { combineGVCFs } from './modules/processGVCFs'
        include { genotypeGVCFs } from './modules/processGVCFs'
        break
    case 'freebayes':
        include { freeBayes } from './modules/freeBayes'
        break
    default: 
        error "Unsupported variant caller: ${params.variant_caller}. Please specify 'haplotype-caller' or 'freebayes'."
}

// Workflow 
workflow {

    // Create qsrc_vcf_ch channel
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

    // Gather read pairs from samplesheet
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) tuple(row[0], [row[1], row[2]])
            else if (row.size() == 3) tuple(row[0], [row[1]])
            else error "Unexpected row format in samplesheet: $row"
        }
    read_pairs_ch.view()

    // Create a channel for the reference genome and prepare .fai and .dict for reference
    reference_ch = Channel.fromPath(params.genome_file)
    prepared_reference_ch = prepareReference(reference_ch)

    // INDEX CHECK & CREATION (FOR ALIGNMENT) 
    // Define paths for index files based on aligner 
    def bwa_index_dir = params.bwa_index_dir
    def bowtie2_index_dir = params.bowtie2_index_dir

    switch(params.aligner) {
        case ['bwa-mem', 'bwa-aln']:
            // Check for existing BWA index files
            def bwa_index_files = new File(bwa_index_dir).listFiles()?.findAll { it.name.endsWith('.amb') }
             if (new File(bwa_index_dir).exists() && bwa_index_files && bwa_index_files.size() > 0) {
                log.info "BWA index files found. Using existing index."
                indexed_genome_ch = Channel.fromPath("${bwa_index_dir}/*")
            } else {
                log.info "No BWA index found. Running bwaIndexGenome module..."
                genome_ch = Channel.fromPath(params.genome_file)
                indexed_genome_ch = bwaIndexGenome(genome_ch)
            }
            break
        case 'bowtie2':
            // Check for existing Bowtie2 index files 
            def bowtie2_index_files = new File(bowtie2_index_dir).listFiles()?.findAll { it.name.endsWith('.bt2') }
            if (new File(bowtie2_index_dir).exists() && bowtie2_index_files && bowtie2_index_files.size() > 0) {
                log.info "Bowtie2 index files found. Using existing index."
                indexed_genome_ch = Channel.fromPath("${bowtie2_index_dir}/*")
            } else {
                log.info "No Bowtie2 index found. Running bowtie2IndexGenome module..."
                genome_ch = Channel.fromPath(params.genome_file)
                indexed_genome_ch = bowtie2IndexGenome(genome_ch)
            }
            break
        default:
            error "Unsupported aligner: ${params.aligner}. Please specify 'bwa-mem', 'bwa-aln', or 'bowtie2'."
    }

    // Run FASTQC on read pairs before fastp
    if (params.fastqc) {
    raw_fastqc_input = read_pairs_ch.map { id, reads ->
        tuple(id, "beforefastp", reads)
    }

    FASTQC_RAW(raw_fastqc_input)
    }

    // Run fastp
    if (params.fastp) {
        trimmed_reads = fastp(read_pairs_ch)
    }

    // Run FASTQC on read pairs after fastp
    if (params.fastqc && params.fastp) {
    trimmed_fastqc_input = trimmed_reads.map { id, reads ->
        tuple(id, "afterfastp", reads)
    }

    FASTQC_TRIMMED(trimmed_fastqc_input)
    }

    // Use the fastp results for alignemnt if fastp is enabled 
    if (params.fastp) {
        read_pairs_for_alignment = fastp.out.fastp_fastqs

        read_pairs_for_alignment.view { id, reads ->
            log.info "ALIGNING (TRIMMED) -> ${id} | ${reads}"
        }
    } else {
        read_pairs_for_alignment = read_pairs_ch

        read_pairs_for_alignment.view { id, reads ->
            log.info "ALIGNING (RAW) -> ${id} | ${reads}"
        }
    }

    // ALIGNMENT 
    def align_ch
    switch(params.aligner) {
        case 'bwa-mem':
            align_ch = alignReadsBwaMem(read_pairs_for_alignment, indexed_genome_ch.collect())
            break
        case 'bwa-aln':
            align_ch = alignReadsBwaAln(read_pairs_for_alignment, indexed_genome_ch.collect())
            break
        case 'bowtie2':
            align_ch = alignReadsBowtie2(read_pairs_for_alignment, indexed_genome_ch.collect())
            break
    }

    // Sort BAM files
    sort_ch = sortBam(align_ch)

    // Mark duplicates in BAM files
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    indexed_bam_ch = indexBam(mark_ch)

    // Conditionally run mapDamage if degraded_dna parameter is set
    if (params.degraded_dna) {
        // Run mapDamage2 process only if degraded_dna is true
        pre_mapDamage_ch = mapDamage2(
            indexed_bam_ch, 
            reference_ch,
        )
        mapDamage_ch = indexMapDamageBam(pre_mapDamage_ch)
    } else {
        // If degraded_dna is not true, just pass through the sorted BAM files
        mapDamage_ch = indexed_bam_ch
    }

    // Create a channel from qsrVcfs
    knownSites_ch = Channel.fromPath(params.qsrVcfs)
        .filter { file -> file.getName().endsWith('.vcf.gz.tbi') || file.getName().endsWith('.vcf.idx') }
        .map { file -> "--known-sites " + file.getBaseName() }
        .collect()

    if (params.bqsr) {
        // Run BQSR on indexed BAM files
        bqsr_ch = baseRecalibrator(
            mapDamage_ch, 
            knownSites_ch, 
            reference_ch, 
            qsrc_vcf_ch.collect()
        )

    } else {
        // If BQSR is skipped, just pass through the mapDamage_ch channel
        bqsr_ch = mapDamage_ch
    }

    // VARIANT CALLING 
    variant_vcf_ch = Channel.empty()

    if (params.variant_caller == "haplotype-caller") {

        gvcf_ch = haplotypeCaller(bqsr_ch, prepared_reference_ch)

        all_gvcf_ch = gvcf_ch
            .toList()
            .map { rows ->
                def sample_ids = rows.collect { it[0] }
                def vcfs = rows.collect { it[1] }
                def idxs = rows.collect { it[2] }
                tuple(sample_ids, vcfs, idxs)
            }

        combined_gvcf_ch = combineGVCFs(all_gvcf_ch, prepared_reference_ch)
        variant_vcf_ch = genotypeGVCFs(combined_gvcf_ch, prepared_reference_ch)

    } else if (params.variant_caller == "freebayes") {

        variant_vcf_ch = freeBayes(bqsr_ch, reference_ch)
    }

    // VARIANT RECALIBRATION OR FILTERING: GATK ONLY 
    filtered_vcf_ch = Channel.empty()

    if (params.variant_recalibration && params.variant_caller == "haplotype-caller") {

        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels': 'known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38': 'known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38': 'known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38': 'known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138': 'known=true,training=false,truth=false,prior=2.0',
            'Mills_and_1000G_gold_standard.indels.hg38': 'known=true,training=true,truth=true,prior=12.0'
        ]

        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { it.getName().endsWith('.vcf.gz') || it.getName().endsWith('.vcf') }
            .map { file ->
                def baseName = file.getName().replaceAll(/\.vcf(\.gz)?$/, '')
                def resourceArgs = resourceOptions.get(baseName) ?: ""
                "--resource:${baseName},${resourceArgs} ${file.getName()}"
            }
            .collect()

        filtered_vcf_ch = variantRecalibrator(
            variant_vcf_ch,
            knownSitesArgs_ch,
            prepared_reference_ch
        )

    } else if (params.variant_caller == "freebayes") {

        // FreeBayes already outputs a final VCF; just pass through
        filtered_vcf_ch = variant_vcf_ch

    } else {

        // fallback: hard filtering OR pass-through
        filtered_vcf_ch = filterVCF(
            variant_vcf_ch,
            prepared_reference_ch
        )
    }

    // Conditionally run identityAnalysis if identity_analysis is true
    if (params.identity_analysis) {

        psam_info_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv(sep: '\t')
            .map { row ->
                if (row.size() == 4) {
                    tuple(row[0], row[3])
                } else if (row.size() == 3) {
                    tuple(row[0], row[2])
                } else {
                    error "Unexpected row format in samplesheet: $row"
                }
            }

        psam_file_ch = psam_info_ch
            .map { sample_info ->
                def sample_id = sample_info[0]
                def sex = sample_info[1] ?: "NA"
                "${sample_id}\t${sample_id}\t0\t0\t${sex}"
            }
            .collect()
            .map { lines ->
                def file = file("${params.outdir}/identity/combined_samples.psam")
                file.parentFile.mkdirs()
                file.text = "#IID\tSID\tPAT\tMAT\tSEX\n" + lines.join("\n") + "\n"
                file
            }

        identity_analysis_ch = identityAnalysis(final_vcf_ch, psam_file_ch)
    }
} 

workflow FASTQC_only {
    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }
    read_pairs_ch.view()

    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}