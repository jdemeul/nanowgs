
/* 
* Polish assembly using Flye
*/
process vcf2diploid {
    label 'cpu_low'
    label 'cpu_low'
    label 'time_low'
    label 'vcf2diploid'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path reads
    path assembly
    val haplotype

    output:
    path "${params.sampleid}_${haplotype}_polished.fasta", emit: polished_assembly

    script:
    """
    flye \
        --polish-target $assembly \
        --nano-hq $reads \
        -t $task.cpus \
        -o ${params.sampleid}_${haplotype}_polished.fasta \
    """

}
