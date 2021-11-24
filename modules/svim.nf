
/* 
* Structural variant calling on aligned reads using SVIM
*/
process svim_sv_calling {
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'svim'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref
    // val step
    // path genomeindex

    output:
    path "working_dir/*.vcf", emit: sv_calls
    path "working_dir/*.png"

    script:
    """
    svim alignment \
        --min_mapq ${params.sv_min_mapq} \
        --min_sv_size ${params.sv_min_size} \
        --sample ${params.sampleid} \
        --read_names \
        working_dir \
        $sorted_bam \
        $genomeref
    """

}
