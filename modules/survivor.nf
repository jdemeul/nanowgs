
/* 
* Creating consensus structural variant calls using SURVIVOR
*/
process survivor_sv_consensus {
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'survivor'

    publishDir path: "${params.outdir}/results/svs_consensus_${step}/", mode: 'copy'

    input:
    path vcfs
    val step

    output:
    path "*_survivor_sv_consensus.vcf", emit: sv_consensus

    script:
    // SURVIVOR merge command structure
    // File with VCF names and paths
    // max distance between breakpoints (0-1 percent of length, 1- number of bp)
    // Minimum number of supporting caller
    // Take the type into account (1==yes, else no)
    // Take the strands of SVs into account (1==yes, else no)
    // Disabled.
    // Minimum size of SVs to be taken into account.
    // Output VCF filename
    """
    ls $vcfs > sv_sample_files
    SURVIVOR merge sv_sample_files \
        ${params.sv_merge_max_dist} \
        2 \
        1 \
        1 \
        0 \
        ${params.sv_min_size} \
        ${params.sampleid}_survivor_sv_consensus.vcf
    SURVIVOR genComp ${params.sampleid}_survivor_sv_consensus.vcf 1 ${params.sampleid}_survivor_sv_consensus.mat.txt
    """
}