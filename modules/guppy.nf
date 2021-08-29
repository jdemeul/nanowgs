
/* 
* Basecall reads with ONT Guppy
*/
process basecall_reads {
    label 'process_middle'
    label 'guppy'
    label ( params.with_gpu ? 'with_gpus': null )

    publishDir path: "${params.outdir}/results/basecalls/", mode: 'copy'

    input:
    path ont_base
    path genomeref
    path genomrefidx

    output:
    path "basecalls/*fastq.gz", emit: fastqs
    path "basecalls/sequencing_summary.txt"

    script:
    """
    guppy_basecaller \
        -i $ont_base \
        -s basecalls \
        -c ${params.guppy_config} \
        --recursive \
        --device "cuda:${params.gpu_devices}" \
        --align_ref $genomeref \
        --compress_fastq \
        --disable_qscore_filtering
    """
}

