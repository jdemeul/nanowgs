
/* 
* Structural variant calling on aligned reads using Sniffles
*/
process sniffles_sv_calling {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'sniffles'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path reference
    // val step

    output:
    path "*_sniffles_svs.vcf.gz", emit: sv_calls
    path "*_sniffles_svs.vcf.gz.tbi", emit: sv_calls_idx
    path "*_sniffles_svs.snf"

    script:
    """
    sniffles --threads $task.cpus \
        --input $sorted_bam \
        --vcf ${params.sampleid}_sniffles_svs.vcf.gz \
        --snf ${params.sampleid}_sniffles_svs.snf \
        --tandem-repeats ${params.tandem_repeats} \
        --minsvlen ${params.sv_min_size} \
        --minsupport ${params.sv_min_support} \
        --reference $reference \
        --output-rnames
        # --phase \
        # --non-germline
        # --mapq ${params.sv_min_mapq}
    """

}
