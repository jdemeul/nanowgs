
/* 
* Assess a de novo genome assembly using Quast
*/
process run_quast {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'quast'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path assembly
    path refgenome
    path refgeneannot
    path ref_bam

    output:
    path "quast_output"

    script:
    """
    quast.py \
        $assembly \
        -r $refgenome \
        -g gene:$refgeneannot \
        -t $task.cpus \
        --large \
        --k-mer-stats \
        --ref-bam $ref_bam \
        -o quast_output
    """
}


   
