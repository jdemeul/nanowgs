
/* 
* Hybrid phasing of SNVs and SVs to create a personal diploid genome reference 
* NOTE UNTESTED IN NEXTFLOW
*/
process create_personal_genome {
    label 'process_low'
    label 'deepvariant'

    publishDir path: "${params.outdir}/results/crossstitch", mode: 'copy'

    input:
    path phased_snps
    path unphased_svs
    path sorted_bam
    path bam_index
    path genomeref
    val karyotype
    val refine

    output:
    path "*.vcf.gz", emit: indel_snv_vcf
    path "*.vcf.gz.tbi", emit: indel_snv_vcf_index
    path "logs"
    path "intermediate_files"
    path "*.visual_report.html"

    script:
    """
    export XDG_CONFIG_HOME="/staging/leuven/stg_00002/lcb/jdemeul/software/"
    module use /staging/leuven/stg_00002/lcb/jdemeul/software/easybuild/modules/all/
    ml --ignore-cache Java/13.0.2 SAMtools/1.9-GCC-6.4.0-2.28 HTSlib/1.9-GCC-6.4.0-2.28 pigz/2.6-GCCcore-6.4.0 Racon/1.4.13-GCCcore-9.3.0

    /staging/leuven/stg_00002/lcb/jdemeul/software/crossstitch/src/crossstitch.sh  \
        /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20.phased_PASS.vcf \
        /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20_LRA_cuteSV_svs_PASS.vcf \
        /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20_minimap2.bam \
        /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chr20.fasta \
        ASA_Edin_BA24_14_18_chr20 \
        xy \
        1
    """
}
