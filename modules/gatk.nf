
/* 
* LiftOver PASS SNPs with GATK
*/
process filter_reads {
    label 'cpu_low'
    label 'mem_low'
    label 'time_mid'
    label 'gatk'

    // publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path genomeref
    path vcf
    path chainfile
    path targetgenomeref

    output:
    path "pass_SNPs_LiftedGRCh28.vcf.gz*"

    script:
    """
    gatk --java-options "-Xmx31G" SelectVariants \
        --exclude-filtered \
        --select-type-to-include SNP \
        -R $genomeref \
        -V $vcf \
        -O pass_SNPs.vcf.gz

    gatk --java-options "-Xmx31G" LiftoverVcf \
        -I pass_SNPs.vcf.gz \
        -O pass_SNPs_LiftedGRCh28.vcf.gz \
        -CHAIN $chainfile \
        -REJECT pass_SNPs_LiftedGRCh28_reject.vcf.gz \
        -R $targetgenomeref \
        --LIFTOVER_MIN_MATCH 0 \
        --RECOVER_SWAPPED_REF_ALT
    """
}
