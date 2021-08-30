
/* 
* SNV and indel calling on aligned reads using Medaka
*/
process medaka_snv_calling {
    label 'process_medium'
    label 'medaka'
    label ("${params.with_gpu}" ? 'with_gpu': null)

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref
    path genomerefidx

    output:
    path "snv_indel_medaka", emit: medaka_variant

    script:
    """
    export CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    medaka_variant \
        -i $sorted_bam \
        -f $genomeref \
        -t $task.cpus \
        -b 150 \
        -n ${params.sampleid} \
        -s ${params.medaka_snp_model} \
        -m ${params.medaka_snp_model.replace("snp", "variant")} \
        -o snv_indel_medaka \
        -p 
    """

}



/* 
* Polish a genome assembly using Medaka
*/
process medaka_assembly_polishing {
    label 'process_high'
    label 'medaka'
    label ("${params.with_gpu}" ? 'with_gpu': null)

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path fastqs
    path draft

    output:
    path "medaka_consensus", emit: consensus

    script:
    """
    export CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    medaka_consensus \
        -i $fastqs \
        -d $draft \
        -o medaka_consensus \
        -m ${medaka_polish_model} \
        -t $task.cpus \
        -b 150
    """

}