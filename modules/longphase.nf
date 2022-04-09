
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
    path "${params.sampleid}_longphase.vcf", emit: snv_indel_phased
    path "${params.sampleid}_longphase_SV.vcf", emit: sv_phased

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
    path "${params.sampleid}_haplotagged.bam", emit: haplotagged_bam
    path hap1ids, emit: hap1ids
    path hap2ids, emit: hap2ids

    shell:
    '''
    longphase haplotag \
        -b !{bam} \
        -t !{task.cpus} \
        -s !{snv_indel} \
        --sv-file !{svs} \
        --percentageThreshold=0.51 \
        --log \
        -o !{params.sampleid}_haplotagged
    
    # generate haplotype files
    grep -v "#" !{params.sampleid}_haplotagged.out | cut -f1,5 | awk '{ $0 = gensub(/\\./, int(rand()*2), "g", $0) }1' > "haplotagged_all.out"
    grep "1$" haplotagged_all.out | cut -f1 > hap1ids
    grep "2$" haplotagged_all.out | cut -f1 > hap2ids
    '''
}