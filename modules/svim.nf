
/* 
* Structural variant calling on aligned reads using SVIM
*/
process svim_sv_calling {
    label 'process_low'
    label 'svim'

    publishDir path: "${params.outdir}/results/svs_svim/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

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
