
/* 
* Assess a de novo genome assembly using Quast
*/
process run_quast {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'quast'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path assembly
    path refgenome
    val haplotype
    // path ref_bam
    // path refgeneannot

    output:
    path "${haplotype}"

    script:
    """
    quast.py \
        $assembly \
        -r $refgenome \
        -t $task.cpus \
        --large \
        --k-mer-stats \
        -o $haplotype
    #    --ref-bam $ref_bam \
    #    -g gene:$refgeneannot \
    """
}


   
