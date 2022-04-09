
/* 
* Sam to sorted bam conversion using samtools
*/
process hapdup {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'hapdup'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path bam
    path bamidx
    path assembly
    // val aligner

    output:
    path "hapdup/*"
    path "hapdup/haplotype_1.fasta", emit: hap1
    path "hapdup/haplotype_2.fasta", emit: hap2

    script:
    """
    hapdup -t $task.cpus \
        --assembly $assembly \
        --bam $bam \
	    --rtype ont \
        --out-dir hapdup
    """

}


/* 
* Generate diploid from a haploid assembly using HapDup while using preassigned haplotagged reads
*/
process hapdup_with_haptagged_bam {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'hapdup'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path haplotaggedbam
    path bamindex
    path assembly

    output:
    path "hapdup/*"
    path "hapdup/haplotype_1.fasta", emit: hap1
    path "hapdup/haplotype_2.fasta", emit: hap2

    script:
    """
    # create files to skip steps1-3 of HapDup pipeline and feed in haplotagged bam
    mkdir -p ./hapdup/margin
    mkdir ./hapdup/pepper
    touch ./hapdup/pepper/PEPPER_VARIANT_FULL.vcf
    cp -P $haplotaggedbam ./hapdup/margin/MARGIN_PHASED.haplotagged.bam
    cp -P $bamindex ./hapdup/margin/MARGIN_PHASED.haplotagged.bam.bai

    hapdup -t $task.cpus \
        --assembly $assembly \
        --bam $haplotaggedbam \
	    --rtype ont \
        --out-dir hapdup
    """

}


/* 
* Polish assembly using Flye
* DEPRECATED AND MOVED TO FLYE MODULE
*/
process flye_polishing {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'hapdup'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path reads
    path assembly
    val haplotype

    output:
    path "${params.sampleid}_${haplotype}_polished.fasta", emit: polished_assembly

    script:
    """
    flye \
        --polish-target $assembly \
        --nano-hq $reads \
        -t $task.cpus \
        -o ${params.sampleid}_${haplotype}_polished.fasta \
    """

}


/* 
* Generate diploid from a haploid assembly using HapDup while using preassigned haplotagged reads
*/
process haptagtransfer {
    label 'cpu_high'
    label 'mem_high'
    label 'time_mid'
    label 'haptagtransfer'

    input:
    path bam
    path assembly

    output:
    path "retagged_sorted.bam", emit: retagged_bam
    path "retagged_sorted.bam.bai", emit: retagged_bamindex

    script:
    """
    echo "Haplotagged bam > haplotagged fastq"
    haplotaggedbam_to_fastq.py -f $bam -o haplotagged_reads.fastq

    echo "Mapping reads to haploid assembly"
    minimap2 -ax map-ont -k 17 -t $task.cpus -L --secondary=no --MD --cap-kalloc=1g -K 10g $assembly haplotagged_reads.fastq > mapped_tagged.sam

    echo "Re-haplotagging reads from read names"
    retag_bam.py -s mapped_tagged.sam -o mapped_tagged.bam

    echo "Sorting retagged bam file"
    samtools sort -@ $task.cpus \
        --write-index \
        -o retagged_sorted.bam##idx##retagged_sorted.bam.bai \
        --reference $assembly \
        -T sorttmp_retagging_sorted \
        mapped_tagged.bam

    """

}