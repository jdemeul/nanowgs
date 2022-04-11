
/* 
* Structural variant calling on aligned reads using cuteSV
*/
process cutesv_sv_calling {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'cutesv'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref
    // val step
    // path genomeindex

    output:
    path "*_cuteSV_svs.vcf", emit: sv_calls

    script:
    """
    cuteSV --min_support ${params.sv_min_support} \
        --report_readid \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --max_cluster_bias_INV 500 \
        --max_cluster_bias_DUP 500 \
        --max_cluster_bias_TRA 50 \
        --diff_ratio_filtering_TRA 0.6 \
        --genotype \
        --min_size ${params.sv_min_size} \
        --sample ${params.sampleid} \
        --threads $task.cpus \
        $sorted_bam \
        $genomeref \
        ${params.sampleid}_cuteSV_svs.vcf \
        `pwd`
        # --min_mapq ${params.sv_min_mapq} \
    """

}