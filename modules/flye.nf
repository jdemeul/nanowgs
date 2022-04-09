
/* 
* Polish assembly using Flye
*/
process flye_polishing {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'flye'

    input:
    path reads
    path assembly
    val haplotype

    output:
    path "polished/polished_${haplotype}.fasta", emit: polished_assembly

    script:
    """
    flye \
        --polish-target $assembly \
        --nano-hq $reads \
        -t $task.cpus \
        -o polished
    mv ./polished/polished_1.fasta ./polished/polished_${haplotype}.fasta
    """

}
