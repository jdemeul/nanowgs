
/* 
* Assess a de novo genome assembly using Quast
*/
process run_mummer {
    label 'cpu_low'
    label 'mem_high'
    label 'time_mid'
    label 'mummer'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path denovoassembly
    path guidedassembly
    val haplotype
    // path ref_bam
    // path refgeneannot

    output:
    path "${haplotype}"

    script:
    """
    dnadiff \
        $guidedassembly \
        $denovoassembly \
        -p $haplotype
    """
}


   
