
/* 
* Sam to sorted bam conversion using samtools
*/
process hapdup {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'hapdup'

    publishDir path: "${params.outdir}/${params.sampleid}/results/${task.process}/", mode: 'copy'

    input:
    path bam
    path bamidx
    path assembly
    // val aligner

    output:
    path "hapdup/*"
    path "hapdup/haplotype_1.fasta", emit: hap1
    path "hapdup/haplotype_2.fasta", emit: hap2

    script:
    """
    hapdup -t $task.cpus \
        --assembly $assembly \
        --bam $bam \
        --out-dir hapdup
    """

}
