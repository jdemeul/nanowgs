
/* 
* Create personal diploid genome using previously phased SNV SV calls 
*/
process create_personal_genome {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'crossstitch'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path phased_snps
    path phased_svs
    path bam
    path genomeref
    val karyotype

    output:
    path "${params.sampleid}.alleleseq/${params.sampleid}.hap1.fa", emit: hap1
    path "${params.sampleid}.alleleseq/${params.sampleid}.hap2.fa", emit: hap2
    path "${params.sampleid}.spliced.scrubbed.vcf.gz"
    path "${params.sampleid}.alleleseq/${params.sampleid}.hap*.chain"

    script:
    """
    /crossstitch/src/crossstitch.sh \
        $phased_snps \
        $phased_svs \
        $bam \
        $genomeref \
        ${params.sampleid} \
        $karyotype \
        0
    """
}


    // /crossstitch/src/crossstitch_short.sh \
    //     $phased_snps \
    //     $phased_svs \
    //     $genomeref \
    //     ${params.sampleid} \
    //     $karyotype




/* 
* Generate diploid from a haploid assembly using HapDup while using preassigned haplotagged reads
*/
process prepare_svs_stitch {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'haptagtransfer'

    input:
    path svs
    path genomeref

    output:
    path "${params.sampleid}_INVaddALT.vcf", emit: fixed_svs

    script:
    """
    echo "Adding ALT allele to inversions"
    fix_inversion_ALT.py --vcf $svs --fasta $genomeref --out ${params.sampleid}_INVaddALT.vcf
    """

}