
/* 
* Structural variant calling on aligned reads using Sniffles
*/
process sniffles_sv_calling {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_low'
    label 'sniffles'

    publishDir path: "${params.outdir}/results/svs_sniffles/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index

    output:
    path "*_sniffles_svs.vcf", emit: sv_calls

    script:
    """
    sniffles --minmapping_qual ${params.sv_min_mapq} \
        --min_support ${params.sv_min_support} \
        --cluster \
        --min_length ${params.sv_min_size} \
        --num_reads_report -1 \
        --threads $task.cpus \
        --tmp_file ./working_dir \
        -m $sorted_bam \
        -v ${params.sampleid}_sniffles_svs.vcf
    """

}
