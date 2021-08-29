
/* 
* Assess a de novo genome assembly using Quast
*/
process run_shasta_assembly {
    label 'process_high'
    label 'quast'

    publishDir path: "${params.outdir}/results/", mode: 'copy'

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


   
