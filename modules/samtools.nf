
/* 
* Sam to sorted bam conversion using samtools
*/
process sam_to_sorted_bam {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'samtools'

    publishDir path: "${params.outdir}/results/bam/", mode: 'copy',
               saveAs: { item -> item.matches("(.*)minimap2(.*)") ? item : null }

    input:
    path mapped_sam
    path genomeref
    // val aligner

    output:
    path "*.bam", emit: sorted_bam
    path "*.bai", emit: bam_index
    path "*stats"

    script:
    def samtools_mem = Math.floor(task.memory.getMega() / task.cpus ) as int
    """
    samtools sort -@ $task.cpus \
        --write-index \
        -o ${params.sampleid}_sorted.bam##idx##${params.sampleid}_sorted.bam.bai \
        -m ${samtools_mem}M \
        --reference $genomeref \
        -T sorttmp_${params.sampleid}_sorted \
        $mapped_sam
    samtools flagstat ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.flagstats
    samtools idxstats ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.idxstats
    # samtools stats ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.stats
    """

}


// /* 
// * Index a fasta file
// */
// process index_fasta {
//     label 'cpu_low'
//     label 'mem_low'
//     label 'time_low'
//     label 'samtools'

//     input:
//     path fasta

//     output:
//     path "*.fai", emit: faidx
//     stdout emit: contig

//     script:
//     """
//     samtools faidx $fasta
//     cat *.fai | cut -f 1
//     """

// }