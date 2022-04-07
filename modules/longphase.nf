
/* 
* Joint phasing of SNVs and SVs using longphase
*/
process longphase_phase {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_low'
    label 'longphase'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path reference
    path snv_indel
    path svs
    path bam
    path bamidx

    output:
    path "${params.sampleid}_longphase*"

    script:
    """
    longphase phase \
        --ont \
        -t $task.cpus \
        -s $snv_indel \
        --sv-file $svs \
        -r $reference \
        -b $bam \
        -o ${params.sampleid}_longphase
    """

}



/* 
* Longphase haplotagging of a bam file
*/
process longphase_tag {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_low'
    label 'longphase'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path snv_indel
    path svs
    path bam
    path bamidx

    output:
    path "${params.sampleid}_haplotagged*"

    script:
    """
    longphase haplotag \
        -b $bam \
        -t $task.cpus \
        -s $snv_indel \
        --sv-file $svs \
        --log \
        -o ${params.sampleid}_haplotagged
    """

}