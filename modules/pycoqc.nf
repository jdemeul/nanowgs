
/* 
* QC a Nanopore run
*/
process run_pycoqc {
    label 'cpu_mid'
    label 'mem_low'
    label 'time_low'
    label 'pycoqc'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path sequencing_summary
    path mapped_bam

    output:
    path "*_pycoQC.html"
    path "*_pycoQC.json"

    script:
    """
    pycoQC -f $sequencing_summary \
        -a $mapped_bam \
        --min_pass_qual ${params.min_read_qscore} \
        --report_title ${params.sampleid} \
        --html_outfile ${params.sampleid}_pycoQC.html \
        --json_outfile ${params.sampleid}_pycoQC.json 
    """
}
