
/* 
* Filtering of SVIM SV calls
*/
process svim_sv_filtering {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path svs
    // val step

    output:
    path "*_Q10.vcf", emit: sv_calls_q10

    script:
    """
    bcftools view -i 'QUAL >= 10' -o ${svs.getSimpleName()}_Q10.vcf $svs
    """
}


/* 
* Filtering of SVIM SV calls
*/
process sniffles_sv_filtering {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path svs
    // val step

    output:
    path "*_PASS.vcf", emit: variants_pass
    path "*.vchk", emit: vcfstats

    script:
    """
    bcftools view -f "PASS" -o ${svs.getSimpleName()}_PASS.vcf $svs
    bcftools stats ${svs.getSimpleName()}_PASS.vcf > ./${svs.getSimpleName()}_PASS.vchk
    """
}


/* 
* Variant call filtering for PASS variants
*/
process variant_filtering {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'
    // ,
        // saveAs: { item -> 
        //                 if ( item.matches("(.*)lra(.*)" ) {
        //                     "/svs/" + item
        //                 } else if item.matches("(.*)medaka(.*)" ) {
        //                     "/snvs_indel_medaka/" + item
        //                 } else {
        //                     "/snvs_indel/" + item
        //                 }}

    input:
    path variants

    output:
    path "${variants.getSimpleName()}_PASS.vcf", emit: variants_pass
    path "${variants.getSimpleName()}_PASS.vchk", emit: vcfstats

    script:
    """
    bcftools view -f "PASS" -e 'ILEN >= 30 | ILEN <= -30' --trim-alt-alleles -o ${variants.getSimpleName()}_PASS.vcf $variants
    bcftools stats ${variants.getSimpleName()}_PASS.vcf > ./${variants.getSimpleName()}_PASS.vchk
    """
}


/* 
* Variant call filtering for PASS variants
*/
process vcf_concat {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path variants

    output:
    path "*_merged_phased_deepvariant.vcf.gz", emit: merged_vcf

    script:
    """
    bcftools concat -O z -o ${params.sampleid}_merged_phased_deepvariant.vcf.gz $variants
    """
}


/* 
* Variant stats
* DEPRECATED incorporated with pass filtering
*/
process vcf_stats {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path variants

    output:
    path "*.vchk", emit: vcfstats

    script:
    """
    bcftools stats -f PASS $variants > ./${variants.getSimpleName()}.vchk
    """
}


/* 
* SV and SNV/indel merging
*/
process vcf_concat_sv_snv {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    input:
    path snv_indels
    path svs

    output:
    path "snv_indel_sv_concat.vcf", emit: merged_vcf

    script:
    """
    echo "${params.sampleid}" > samples
    bcftools reheader -s samples $snv_indels | bcftools view -o reheadered_snv_indels.vcf.gz
    bcftools reheader -s samples $svs | bcftools view -o reheadered_svs.vcf.gz 
    bcftools index reheadered_snv_indels.vcf.gz
    bcftools index reheadered_svs.vcf.gz
    bcftools concat -a -o snv_indel_sv_concat.vcf reheadered_snv_indels.vcf.gz reheadered_svs.vcf.gz
    """
}
