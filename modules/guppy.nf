
/* 
* Basecall reads with ONT Guppy
*/
process basecall_reads {
    label 'guppy'
    label ( params.with_gpu ? 'with_gpu': 'cpu_high')
    label ( params.with_gpu ? null: 'mem_mid')
    label ( params.with_gpu ? null: 'time_high')

    publishDir path: "${params.outdir}/results/basecalls/", mode: 'copy'

    input:
    path ont_base
    path genomeref
    // path genomrefidx

    output:
    path "basecalls/*fastq.gz", emit: fastqs
    path "basecalls/sequencing_summary.txt"

    script:
    if ( params.with_gpu )
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
    else 
        """
        guppy_basecaller \
            -i $ont_base \
            -s basecalls \
            --recursive \
            --device "cuda:${params.gpu_devices}" \
            --align_ref $genomeref \
            --compress_fastq \
            --disable_qscore_filtering
        """
}

