
/* 
* SNV and indel calling on aligned reads using Medaka
*/
process medaka_snv_calling {
    label 'medaka'
    label ( params.with_gpu ? 'with_gpu': 'cpu_high')
    label ( params.with_gpu ? null: 'mem_high')
    label ( params.with_gpu ? null: 'time_high')

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref
    // path genomerefidx

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
    label 'medaka'
    label ( params.with_gpu ? 'with_gpu': 'cpu_low')
    label ( params.with_gpu ? null: 'mem_mid')
    label ( params.with_gpu ? null: 'time_high')

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path fastqs
    path draft

    output:
    path "medaka_consensus/consensus.fasta", emit: consensus

    script:
    """
    export CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    medaka_consensus \
        -i $fastqs \
        -d $draft \
        -o medaka_consensus \
        -m ${params.medaka_polish_model} \
        -t $task.cpus \
        -b 150
    """

}



/* 
* Polish a genome assembly using Medaka
*/
process medaka_assembly_polish_align {
    label 'medaka'
    label 'cpu_high'
    label 'mem_mid'
    label 'time_mid'

    // publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path fastqs
    path draft

    output:
    path "calls_to_draft.bam", emit: calls_to_draft
    path "calls_to_draft.bam.bai", emit: calls_to_draft_index

    script:
    """
    mini_align \
        -i $fastqs \
        -r $draft \
        -m \
        -p calls_to_draft \
        -t $task.cpus
    """

}


/* 
* Polish a genome assembly using Medaka
*/
process medaka_assembly_polish_consensus {
    label 'medaka'
    label 'cpu_low'
    // cpus = 4
    label 'mem_mid'
    label 'time_low'

    input:
    path bam
    path bamidx
    each contig

    output:
    path "*.hdf", emit: probs

    script:
    """
    # export CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    medaka consensus \
        $bam \
        --model ${params.medaka_polish_model} \
        ${contig[0]}.hdf \
        --batch 200 \
        --threads $task.cpus \
        --regions ${contig.join(" ")}
    """

}


/* 
* Polish a genome assembly using Medaka
*/
process medaka_assembly_polish_stitch {
    label 'medaka'
    label 'cpu_low'
    label 'mem_low'
    label 'time_low'

    publishDir path: "${params.outdir}/results/racon_medaka_consensus/", mode: 'copy'

    input:
    path hdfs
    path draft

    output:
    path "consensus.fasta", emit: consensus

    script:
    """
    medaka stitch \
        --threads $task.cpus \
        $hdfs \
        $draft \
        consensus.fasta 
    """

}