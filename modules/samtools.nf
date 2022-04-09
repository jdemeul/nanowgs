
/* 
* Sam to sorted bam conversion using samtools
*/
process sam_to_sorted_bam {
    label 'cpu_mid'
    label 'mem_mid'
    label 'time_mid'
    label 'samtools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'
    // publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy',
    //            saveAs: { item -> item.matches("(.*)stats") ? item : null }

    input:
    path mapped_sam
    path genomeref
    // val aligner

    output:
    path "${params.sampleid}_sorted.bam", emit: sorted_bam
    path "${params.sampleid}_sorted.bam.bai", emit: bam_index
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
    samtools idxstats -@ 4 ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.idxstats
    samtools stats -@ 4 ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.stats
    """

}



/* 
* Sam to sorted bam conversion using samtools
*/
process get_haplotype_readids {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'samtools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path haplotagged_bam

    output:
    path "hap1ids", emit: hap1ids
    path "hap2ids", emit: hap2ids
    path "hap0ids"

    script:
    """
    samtools view -@ $task.cpus -d "HP:0" $haplotagged_bam | cut -f 1 | shuf > hap0ids
    split -n l/2 hap0ids
    mv xaa hap1ids
    mv xab hap2ids
    samtools view -@ $task.cpus -d "HP:1" $haplotagged_bam | cut -f 1 >> hap1ids
    samtools view -@ $task.cpus -d "HP:2" $haplotagged_bam | cut -f 1 >> hap2ids
    """

}


/* 
* Sam to sorted bam conversion using samtools
*/
process index_bam {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'samtools'

    input:
    path mapped_bam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bam_index

    """
    samtools index -b -@ $task.cpus $mapped_bam
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