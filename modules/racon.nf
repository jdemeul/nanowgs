
/* 
* Polish a genome assembly using Racon
*/
process racon_assembly_polishing {
    label 'racon'
    // label ("${params.with_gpu}" ? 'with_gpu': null)
    label 'bigmem'

    // publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path fastq
    path aligned_reads_sam
    path draft

    output:
    path "consensus.fa", emit: consensus 
    // stdout emit: consensus

    script:
    """
    # export CUDA_VISIBLE_DEVICES=${params.gpu_devices} memory requirements do not permit running on the GPU nodes
    racon \
    -u \
    -m 8 -x -6 -g -8 -w 500 \
    -t 28 \
    $fastq $aligned_reads_sam $draft > consensus.fa
    echo "Done running Racon and output to consensus.fa"
    # --cudapoa-batches 50 \  
    # --cudaaligner-batches 10 \
    """

    stub:
    """
    echo "racon -u -m 8 -x -6 -g -8 -w 500 -t $task.numthreads $fastq $aligned_reads_sam $draft > racon_consensus.fa"
    touch racon_consensus.fa
    """

}