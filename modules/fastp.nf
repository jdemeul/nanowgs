
/* 
* Process (filter + trim) fastq files with fastp
*/
process filter_reads {
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'fastp'

    publishDir path: "${params.outdir}/results/fastq/", mode: 'copy'

    input:
    path fastqs

    output:
    path "*_trimmed.fastq.gz", emit: fastq_trimmed
    path "fastp*"

    script:
    """
    zcat $fastqs | fastp --stdin \
        --disable_adapter_trimming \
        --average_qual ${params.min_read_qscore} --qualified_quality_phred 0 --unqualified_percent_limit 100 --n_base_limit 50 \
        --length_required 100 \
        --trim_front1 30 --trim_tail1 15 \
        -o ${params.sampleid}_trimmed.fastq.gz \
        --thread $task.cpus 
    """
}
