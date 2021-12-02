
/* 
* Sam to sorted bam conversion using samtools
*/
process dipdiff {
    label 'cpu_high'
    label 'mem_high'
    label 'time_low'
    label 'dipdiff'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path reference
    path hap1_fasta
    path hap2_fasta

    output:
    path "dipdiff"

    script:
    """
    dipdiff.py -t 18 \
        --reference $reference \
        --pat $hap1_fasta \
        --mat $hap2_fasta \
        --out-dir dipdiff
    """

}
