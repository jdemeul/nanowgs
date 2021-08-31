
/* 
* Polish a genome assembly using Racon
*/
process racon_assembly_polishing {
    label 'process_high'
    label 'racon'
    label ("${params.with_gpu}" ? 'with_gpu': null)

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path fastq
    path aligned_reads_sam
    path draft

    output:
    path "consensus.fasta", emit: consensus

    script:
    """
    export CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    racon \
    -t $task.cpus \
    --cudapoa-batches 50 \
    --cudaaligner-batches 10 \
    $fastq \
    $aligned_reads_sam \
    $draft
    """

}