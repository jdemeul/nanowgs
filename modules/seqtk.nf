
/* 
* Joint phasing of SNVs and SVs using longphase
*/
process seqtk {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'seqtk'

    input:
    path hap1ids
    path hap2ids
    path reads

    output:
    path "hap1reads.fq", emit: hap1reads
    path "hap2reads.fq", emit: hap2reads

    script:
    """
    echo "Collecting HP1 reads"
    seqtk subseq $reads $hap1ids > hap1reads.fq
    echo "Collecting HP2 reads"
    seqtk subseq $reads $hap2ids > hap2reads.fq
    """

}
