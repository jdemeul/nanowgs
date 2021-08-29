
/* 
* Structural variant calling on aligned reads using cuteSV
*/
process cutesv_sv_calling {
    label 'process_high'
    label 'cutesv'

    publishDir path: "${params.outdir}/results/svs_cutesv/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "*_cuteSV_svs.vcf", emit: sv_calls

    script:
    """
    cuteSV --min_support ${params.sv_min_support} \
        --min_mapq ${params.sv_min_mapq} \
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
    """

}