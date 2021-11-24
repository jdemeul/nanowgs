
/* 
* Basecall reads with ONT Guppy
*/
process basecall_reads {
    label 'guppy'
    label 'with_p100node'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path ont_base
    path genomeref
    // path genomrefidx

    output:
    path "basecalls/*fastq.gz", emit: fastqs
    path "basecalls/sequencing_summary.txt"

    script:
    """
    /opt/ont-guppy/guppy_basecaller \
        -i $ont_base \
        -s ./basecalls \
        -c ${params.guppy_config} \
        --recursive \
        --device "cuda:all" \
        --align_ref $genomeref \
        --compress_fastq \
        --disable_qscore_filtering
    """
}

