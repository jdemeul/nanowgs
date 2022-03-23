
/* 
* Structural variant calling on aligned reads using Dysgu
*/
process dysgu_sv_calling {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'dysgu'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path reference
    // val step

    output:
    path "*_dysgu_svs.vcf", emit: sv_calls

    script:
    """
    dysgu run \
        --mode nanopore \
        --diploid True \
        --min-support ${params.sv_min_support} \
        --min-size ${params.sv_min_size} \
        --max-cov auto \
        --mq ${params.sv_min_mapq} \
        -o ./${params.sampleid}_dysgu_svs.vcf \
        -p $task.cpus \
        -c \
        $reference \
        ./dysgu/ \
        $sorted_bam
    """

}
