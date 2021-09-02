
/* 
* Filtering of SVIM SV calls
*/
process svim_sv_filtering {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'bcftools'

    publishDir path: "${params.outdir}/results/svs_svim/", mode: 'copy'

    input:
    path svs

    output:
    path "*_Q10.vcf", emit: sv_calls_q10

    script:
    """
    bcftools view -i 'QUAL >= 10' -o ${svs.getSimpleName()}_Q10.vcf $svs
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

    publishDir path: "${params.outdir}/results/", mode: 'copy'
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
    path "*_PASS.vcf", emit: variants_pass

    script:
    """
    bcftools view -f "PASS" -o ${variants.getSimpleName()}_PASS.vcf $variants
    """
}