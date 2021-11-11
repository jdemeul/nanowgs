
/* 
* Sam to sorted bam conversion using samtools
*/
process dipdiff {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'dipdiff'

    publishDir path: "${params.outdir}/${params.sampleid}/results/${task.process}/", mode: 'copy'

    input:
    path reference
    path hap1_fasta
    path hap2_fasta

    output:
    path ""

    script:
    """
    dipdiff.py -t $task.cpus \
        --reference $reference \
        --pat $hap1_fasta \
        --mat $hap2_fasta \
        --out-dir dipdiff
    """

}
