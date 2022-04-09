
/* 
* Process (filter + trim) fastq files with fastp
*/
process parallel_gzip {
    label 'cpu_mid'
    label 'mem_low'
    label 'time_low'
    label 'pigz'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path fastq

    output:
    path "${params.sampleid}.gz", emit: fastqgz

    script:
    """
    pigz \
        --best \
        --stdout \
        --processes $task.cpus \
        $fastq > ${params.sampleid}.gz
    """
}
