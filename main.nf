#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

      nextflow run khigashi1987/CUTRUN_Nextflow --input ./samples.csv \
        --adapter /path/to/adapter.fa \
        --fasta /path/to/genome.fa \
        --gtf /path/to/genes.gtf \
        --bt2_index /path/to/genome \
        --gene_bed /path/to/genes.bed \
        --macs_gsize "2.7e9" \
        -profile docker

    Mandatory arguments:
      --input [file]                  Comma-separated file containing information about the samples in the experiment
      --adapter [file]                Path to Fasta file of adapter sequences to be trimmed.
      --fasta [file]                  Path to Fasta reference genome sequence.
      --gtf [file]                    Path to GTF file.
      --bt2_index [file]              Full path to directory containing Bowtie2 index including base name i.e. /path/to/index/genome
      --gene_bed [file]               Path to BED file containing gene intervals
      -profile [str]                  Configuration profile to use. Available: docker

    Preprocessing
      --blacklist [file]              Path to blacklist regions (.BED format), used for filtering peaks

    Trimming
      --seq_len [int]                 length of fastq sequences. (Default: 100)

    Alignment
      --skip_dedup [bool]             Skip duplicate removal process. (Default: false)
      --skip_filter [bool]            Skip filtering process of large fragments. (Default: false)

    Peaks
      --macs_gsize [str]              Effective genome size of reference. (e.g., "2.7e9" for hg38, "1.87e9" for mm10)
      --macs_qvalue [float]           Minimum FDR (q-value) cutoff for peak detection. (Default: 0.01)
      --skip_seacr [bool]             Skip SEACR-peak calling process. (Default: false)
      --seacr_mode [str]              SEACR running mode. ["relaxed" | "stringent"] (Default: "stringent")
      --seacr_threshold [float]       A numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.
    
    Plotting
      --skip_plot_profile [bool]      Skip PLOTPROFILE process. (Default: false)
    
    QC
      --skip_fastqc [bool]            Skip FastQC. (Default: false)
      --skip_multiqc [bool]           Skip MultiQC. (Default: false)

    Other
      --outdir [file]                 The output directory where the results will be saved (Default: './results')
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic (Default: false)
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

////////////////////////////////////////////////////
/* --        MULTIQC CONFIG FILES              -- */
////////////////////////////////////////////////////

// Pipeline config
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

// Header files for MultiQC
ch_peak_count_header = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

Channel.fromPath(params.input, checkIfExists: true)
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
    .into {  ch_raw_reads_fastqc;
            ch_raw_reads_trim }

if (params.gtf)       { ch_gtf = file(params.gtf, checkIfExists: true) } else { exit 1, 'GTF annotation file not specified!' }
if (params.gene_bed)  { ch_gene_bed = file(params.gene_bed, checkIfExists: true) }
if (params.blacklist) { ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true) } else { ch_blacklist = Channel.empty() }

if (params.fasta) {
    ch_fasta = file(params.fasta, checkIfExists: true)
} else {
    exit 1, 'Fasta file not specified!'
}

if (params.adapter) {
    ch_adapter = file(params.adapter, checkIfExists: true)
} else {
    exit 1, 'Adapter file not specified!'
}

if (params.bt2_index) {
    lastPath = params.bt2_index.lastIndexOf(File.separator)
    bt2_dir  = params.bt2_index.substring(0,lastPath+1)
    bt2_base = params.bt2_index.substring(lastPath+1)
    Channel
        .fromPath(bt2_dir, checkIfExists: true)
        .set { ch_bt2_index }
}

process MAKE_GENOME_FILTER {
    tag "$fasta"
    publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

    input:
    path fasta from ch_fasta
    path blacklist from ch_blacklist.ifEmpty([])

    output:
    path "$fasta"                                      // FASTA FILE
    path '*.fai'                                       // FAI INDEX FOR REFERENCE GENOME
    path '*.bed' into ch_genome_filter_regions         // BED FILE WITHOUT BLACKLIST REGIONS
    path '*.sizes' into ch_genome_sizes_bigwig,        // CHROMOSOME SIZES FILE FOR BEDTOOLS
                        ch_genome_sizes_bedgraph

    script:
    blacklist_filter = params.blacklist ? "sortBed -i $blacklist -g ${fasta}.sizes | complementBed -i stdin -g ${fasta}.sizes" : "awk '{print \$1, '0' , \$2}' OFS='\t' ${fasta}.sizes"
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    $blacklist_filter > ${fasta}.include_regions.bed
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                       HEADER LOG INFO                               -- */
///////////////////////////////////////////////////////////////////////////////

// Header log info
def summary = [:]
summary['Run Name']               = params.name
summary['Input File']             = params.input
summary['Adapter File']           = params.adapter
summary['Fasta File']             = params.fasta
summary['GTF File']               = params.gtf
summary['Gene BED File']          = params.gene_bed
summary['Bowtie2 Index']          = params.bt2_index
summary['Blacklist BED']          = params.blacklist
summary['MACS2 Genome Size']      = params.macs_gsize ?: 'Not supplied.'
if (params.macs_gsize)            summary['MACS2 q-value'] = params.macs_qvalue
if (!params.skip_seacr)           summary['SEACR Mode'] = params.seacr_mode
if (!params.skip_seacr)           summary['SEACR Threshold'] = params.seacr_threshold
if (params.skip_fastqc)           summary['Skip FastQC'] = 'Yes'
if (params.skip_multiqc)          summary['Skip MultiQC'] = 'Yes'
if (workflow.containerEngine)     summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output Dir']             = params.outdir
summary['Launch Dir']             = workflow.launchDir
summary['Working Dir']            = workflow.workDir
summary['Script Dir']             = workflow.projectDir
summary['User']                   = workflow.userName
summary['Config Profile']         = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join('\n')

///////////////////////////////////////////////////////////////////////////////
/* --                        FASTQ QC                                     -- */
///////////////////////////////////////////////////////////////////////////////

process FASTQC {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.endsWith('.zip') ? "zips/$filename" : filename
                }

    when:
    !params.skip_fastqc

    input:
    tuple val(name), path(reads) from ch_raw_reads_fastqc

    output:
    path '*.{zip,html}' into ch_fastqc_reports_mqc

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    """
    [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
    [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
    fastqc -q -t $task.cpus ${name}_1.fastq.gz
    fastqc -q -t $task.cpus ${name}_2.fastq.gz
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                        ADAPTER TRIMMING                             -- */
///////////////////////////////////////////////////////////////////////////////

process TRIM_1STROUND {
    tag "$name"
    publishDir "${params.outdir}/trim_1st", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith('.html')) "fastqc/$filename"
                    else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
                    else filename
            }

    input:
    tuple val(name), path(reads) from ch_raw_reads_trim
    path adapter from ch_adapter

    output:
    tuple val(name), path('*.paired.fastq.gz') into ch_trimmed_1st
    path '*.{zip,html}' into ch_trim1st_fastqc_reports_mqc

    """
    [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
    [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
    trimmomatic PE -threads 1 -phred33 \
      ${name}_1.fastq.gz ${name}_2.fastq.gz \
      ${name}_1.paired.fastq.gz ${name}_1.unpaired.fastq.gz \
      ${name}_2.paired.fastq.gz ${name}_2.unpaired.fastq.gz \
      ILLUMINACLIP:${adapter}:2:15:4:4:true \
      LEADING:20 \
      TRAILING:20 \
      SLIDINGWINDOW:4:15 \
      MINLEN:25
    fastqc -q -t $task.cpus ${name}_1.paired.fastq.gz
    fastqc -q -t $task.cpus ${name}_2.paired.fastq.gz
    """
}

process TRIM_2NDROUND {
    tag "$name"
    publishDir "${params.outdir}/trim_2nd", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith('.html')) "fastqc/$filename"
                    else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
                    else filename
            }

    input:
    tuple val(name), path(reads) from ch_trimmed_1st

    output:
    tuple val(name), path('*.paired.trimmed.fastq.gz') into ch_trimmed_2nd
    path '*.{zip,html}' into ch_trim2nd_fastqc_reports_mqc

    """
    kseq_test ${reads[0]} ${params.seq_len} ${name}_1.paired.trimmed.fastq.gz
    kseq_test ${reads[1]} ${params.seq_len} ${name}_2.paired.trimmed.fastq.gz
    fastqc -q -t $task.cpus ${name}_1.paired.trimmed.fastq.gz
    fastqc -q -t $task.cpus ${name}_2.paired.trimmed.fastq.gz
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                           ALIGN                                     -- */
///////////////////////////////////////////////////////////////////////////////

process BOWTIE2 {
    tag "$name"
    publishDir "${params.outdir}/bowtie2_aln", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads) from ch_trimmed_2nd
    path index from ch_bt2_index.collect()

    output:
    tuple val(name), path('*.bam') into ch_bt2

    """
    bowtie2 -p 2 --dovetail --phred33 -x ${index}/${bt2_base} \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
    | samtools view -bS - -o ${name}.bam \
    """
}

process SORT_BAM {
    tag "$name"
    if (params.save_align_intermeds) {
        publishDir "${params.outdir}/sorted_bam", mode: params.publish_dir_mode
    }

    input:
    tuple val(name), path(bam) from ch_bt2

    output:
    tuple val(name), path('*.sorted.bam') into ch_sorted_bam

    """
    samtools view -bh -f 3 -F 4 -F 8 ${bam} > ${name}.filter_unmapped.bam
    java -jar ${params.picard_jar} SortSam INPUT=${name}.filter_unmapped.bam \
      OUTPUT=${name}.sorted.bam \
      SORT_ORDER=coordinate \
      VALIDATION_STRINGENCY=SILENT
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                      FILTER  ALIGNMENT                              -- */
///////////////////////////////////////////////////////////////////////////////

if (params.skip_dedup) {
    ch_dedup_bam = ch_sorted_bam
} else {
    process DEDUP {
      tag "$name"
      publishDir "${params.outdir}/remove_duplicates", mode: params.publish_dir_mode

      input:
      tuple val(name), path(bam) from ch_sorted_bam

      output:
      tuple val(name), path('*.sorted.dedup.bam') into ch_dedup_bam
      path '*.txt' into ch_picard_metrics_mqc

      """
      java -jar ${params.picard_jar} MarkDuplicates INPUT=${bam} \
        OUTPUT=${name}.sorted.marked.bam \
        VALIDATION_STRINGENCY=SILENT \
        METRICS_FILE=metrics.${name}.txt
      samtools view -bh -F 1024 ${name}.sorted.marked.bam > ${name}.sorted.dedup.bam
      """
    }
}

if (params.skip_filter) {
    ch_filtered_bam = ch_dedup_bam
} else {
    process FILTER120 {
        tag "$name"
        publishDir "${params.outdir}/filter120", mode: params.publish_dir_mode

        input:
        tuple val(name), path(bam) from ch_dedup_bam

        output:
        tuple val(name), path('*.sorted.dedup.filtered.{bam,bam.bai}') into ch_filtered_bam

        """
        samtools view -h ${bam} \
          |LC_ALL=C awk -f /cutruntools/filter_below.awk \
          |samtools view -Sb - -o ${name}.sorted.dedup.filtered.bam
        samtools index ${name}.sorted.dedup.filtered.bam
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
/* --                      CONVERT  ALIGNMENT                             -- */
///////////////////////////////////////////////////////////////////////////////

ch_filtered_bam.into {ch_filtered_bam_macs; ch_filtered_bam_filter_macs; ch_filtered_bam_flagstat; ch_filtered_bam_bedgraph; ch_filtered_bam_bigwig}

process FLAGSTAT {
    tag "$name"
    publishDir "${params.outdir}/flagstat", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam) from ch_filtered_bam_flagstat

    output:
    tuple val(name), path('*.flagstat') into ch_flagstat_bigwig,
                                             ch_flagstat_macs
    path '*.{flagstat,idxstats,stats}' into ch_flagstat_mqc

    """
    samtools flagstat ${bam[0]} > ${name}.sorted.dedup.filtered.bam.flagstat
    samtools idxstats ${bam[0]} > ${name}.sorted.dedup.filtered.bam.idxstats
    samtools stats ${bam[0]} > ${name}.sorted.dedup.filtered.bam.stats
    """
}

// Read depth normalised bigWig
process BIGWIG {
    tag "$name"
    publishDir "${params.outdir}/bigwig", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam), path(flagstat) from ch_filtered_bam_bigwig.join(ch_flagstat_bigwig, by: [0])
    path sizes from ch_genome_sizes_bigwig.collect()

    output:
    tuple val(name), path('*.bigWig') into ch_bigwig_plotprofile
    path '*scale_factor.txt'

    """
    SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${name}.scale_factor.txt
    genomeCoverageBed -ibam ${bam[0]} -bg -scale \$SCALE_FACTOR -pc | sort -T '.' -k1,1 -k2,2n > ${name}.bedGraph
    bedGraphToBigWig ${name}.bedGraph $sizes ${name}.bigWig
    find * -type f -name "*.bigWig" -exec echo -e "bigwig/"{}"\\t0,0,178" \\; > ${name}.bigWig.igv.txt
    """
}

// BAM to Bedgraph
process BEDGRAPH {
    tag "$name"
    publishDir "${params.outdir}/bedgraph", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam) from ch_filtered_bam_bedgraph
    path sizes from ch_genome_sizes_bedgraph.collect()

    output:
    tuple val(name), path('*.bedgraph') into ch_bedgraph_seacr
    path '*.bedgraph'

    """
    bedtools bamtobed -bedpe -i ${bam[0]} > ${name}.bed
    awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${name}.bed > ${name}.clean.bed
    cut -f 1,2,6 ${name}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${name}.fragments.bed
    bedtools genomecov -bg -i ${name}.fragments.bed -g ${sizes} > ${name}.fragments.bedgraph
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                          PEAK CALL                                  -- */
///////////////////////////////////////////////////////////////////////////////

process MACS2 {
    tag "$name"
    publishDir "${params.outdir}/macs", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith('.tsv')) "qc/$filename"
                      else if (filename.endsWith('.igv.txt')) null
                      else filename
                }

    when:
    params.macs_gsize

    input:
    tuple val(name), path(bam) from ch_filtered_bam_macs

    output:
    tuple val(name), path('*.narrowPeak') into ch_macs_peaks
    path '*.{bed,xls,bdg}'

    """
    macs2 callpeak \
      -t ${bam[0]} \
      -g $params.macs_gsize \
      -f BAMPE \
      -n ${name} \
      -q ${params.macs_qvalue} \
      -B --SPMR \
      --keep-dup all 
    """
}

process SEACR {
    tag "$name"
    publishDir "${params.outdir}/seacr", mode: params.publish_dir_mode

    when:
    !params.skip_seacr

    input:
    tuple val(name), path(bdg) from ch_bedgraph_seacr

    output:
    tuple val(name), path('*.narrowPeak') into ch_seacr_peaks
    path '*.narrowPeak'

    """
    bash /SEACR/SEACR_1.3.sh \
      ${bdg} \
      ${params.seacr_threshold} \
      non \
      ${params.seacr_mode} \
      ${name}
    awk 'BEGIN {FS="\t"; OFS="\t";} {print \$1,\$2,\$3,\$6,\$4,".",\$5}' ${name}.${params.seacr_mode}.bed \
      > ${name}.${params.seacr_mode}.narrowPeak
    """
}

if (params.macs_gsize && !params.skip_seacr) {
    ch_genome_filter_regions.into {ch_genome_filter_regions_macs; ch_genome_filter_regions_seacr}
} else if (params.macs_gsize) {
    ch_genome_filter_regions_macs = ch_genome_filter_regions
} else {
    ch_genome_filter_regions_seacr = ch_genome_filter_regions
}

if (params.macs_gsize && !params.blacklist) {
    ch_macs_peaks.into {ch_filtered_macs_peaks_qc; ch_filtered_macs_peaks_annot}
} else {
    process FILTER_MACS_PEAKS {
        tag "$name"
        publishDir "${params.outdir}/macs", mode: params.publish_dir_mode

        input:
        tuple val(name), path(peaks), path(flagstat), path(bam) from ch_macs_peaks.join(ch_flagstat_macs, by: [0]).join(ch_filtered_bam_filter_macs, by:[0])
        path(bed) from ch_genome_filter_regions_macs.collect()
        path peak_count_header from ch_peak_count_header
        path frip_score_header from ch_frip_score_header


        output:
        tuple val(name), path('*.filtered.narrowPeak') into ch_filtered_macs_peaks_qc,
                                                            ch_filtered_macs_peaks_annot
        path '*_mqc.tsv' into ch_macs_mqc

        """
        cat ${peaks} | \
          grep -v -e "chrM" | \
          sort-bed - | \
          bedops -n 1 - ${bed} > ${name}.filtered.narrowPeak

        cat ${name}.filtered.narrowPeak | \
            wc -l | \
            awk -v OFS='\t' '{ print "${name}", \$1 }' | \
            cat $peak_count_header - > ${name}_peaks.count_mqc.tsv
        READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${name}.filtered.narrowPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
        grep 'mapped (' $flagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${name}", a/\$1}' | \
            cat $frip_score_header - \
            > ${name}_peaks.FRiP_mqc.tsv
        """
    }
}

if (!params.skip_seacr && !params.blacklist) {
    ch_seacr_peaks.into {ch_filtered_seacr_peaks; ch_filtered_seacr_peaks_annot}
} else {
    process FILTER_SEACR_PEAKS {
        tag "$name"
        publishDir "${params.outdir}/seacr", mode: params.publish_dir_mode

        input:
        tuple val(name), path(peaks) from ch_seacr_peaks
        path(bed) from ch_genome_filter_regions_seacr.collect()

        output:
        tuple val(name), path('*.filtered.narrowPeak') into ch_filtered_seacr_peaks,
                                                     ch_filtered_seacr_peaks_annot

        """
        cat ${peaks} | \
          grep -v -e "chrM" | \
          sort-bed - | \
          bedops -n 1 - ${bed} > ${name}.${params.seacr_mode}.filtered.narrowPeak
        """
    }
}

process PEAK_ANNOTATION_MACS {
    tag "$name"
    publishDir "${params.outdir}/macs", mode: params.publish_dir_mode

    when:
    params.macs_gsize

    input:
    tuple val(name), path(peaks) from ch_filtered_macs_peaks_annot
    path fasta from ch_fasta
    path gtf from ch_gtf

    output:
    path '*.annotatePeaks.txt'

    """
    annotatePeaks.pl \
        ${peaks} \
        ${fasta} \
        -gid \
        -gtf ${gtf} \
        > ${name}.annotatePeaks.txt
    """
}

process PEAK_ANNOTATION_SEACR {
    tag "$name"
    publishDir "${params.outdir}/seacr", mode: params.publish_dir_mode

    when:
    !params.skip_seacr

    input:
    tuple val(name), path(peaks) from ch_filtered_seacr_peaks_annot
    path fasta from ch_fasta
    path gtf from ch_gtf

    output:
    path '*.annotatePeaks.txt'

    """
    annotatePeaks.pl \
        ${peaks} \
        ${fasta} \
        -gid \
        -gtf ${gtf} \
        > ${name}.annotatePeaks.txt
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                          PLOTTING                                   -- */
///////////////////////////////////////////////////////////////////////////////

process PLOTPROFILE {
    tag "$name"
    publishDir "${params.outdir}/deepTools/plotProfile", mode: params.publish_dir_mode

    when:
    !params.skip_plot_profile

    input:
    tuple val(name), path(bigwig) from ch_bigwig_plotprofile
    path bed from ch_gene_bed

    output:
    path '*.plotProfile.tab' into ch_plotprofile_mqc
    path '*.{gz,pdf,mat.tab}'

    script:
    """
    computeMatrix scale-regions \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${name}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${name}.computeMatrix.vals.mat.tab \\
        --regionBodyLength 1000 \\
        --beforeRegionStartLength 3000 \\
        --afterRegionStartLength 3000 \\
        --skipZeros \\
        --smartLabels \\
        --numberOfProcessors $task.cpus
    plotProfile --matrixFile ${name}.computeMatrix.mat.gz \\
        --outFileName ${name}.plotProfile.pdf \\
        --outFileNameData ${name}.plotProfile.tab
    plotHeatmap --matrixFile ${name}.computeMatrix.mat.gz \\
        --outFileName ${name}.plotHeatmap.pdf \\
        --outFileNameMatrix ${name}.plotHeatmap.mat.tab
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                          MULTIQC                                    -- */
///////////////////////////////////////////////////////////////////////////////

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'cutrun-nextflow-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'cutrun_nextflow Workflow Summary'
    section_href: 'https://github.com/khigashi1987/CUTRUN_Nextflow'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

process MULTIQC {
    tag "$name"
    publishDir "${params.outdir}/multiqc/", mode: params.publish_dir_mode

    when:
    !params.skip_multiqc && params.macs_gsize && params.blacklist

    input:
    path (multiqc_config) from ch_multiqc_config
    path workflow_summary from ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    path ('fastqc/*') from ch_fastqc_reports_mqc.collect().ifEmpty([])
    path ('trim_1st/fastqc/*') from ch_trim1st_fastqc_reports_mqc.collect().ifEmpty([])
    path ('trim_2nd/fastqc/*') from ch_trim2nd_fastqc_reports_mqc.collect().ifEmpty([])
    path ('alignment/*') from ch_flagstat_mqc.collect()
    path ('alignment/picard_metrics/*') from ch_picard_metrics_mqc.collect()
    path ('macs/*') from ch_macs_mqc.collect().ifEmpty([])
    path ('deeptools/*') from ch_plotprofile_mqc.collect().ifEmpty([])

    output:
    path '*multiqc_report.html'
    path '*_data'

    script:
    rtitle = "--title \"$params.name\""
    rfilename = "--filename " + params.name + "_multiqc_report"
    """
    multiqc . -f $rtitle $rfilename
    """
}