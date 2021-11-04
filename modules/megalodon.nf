
/* 
* Identification of modified bases using Megalodon
*/
process megalodon {
    label 'megalodon'
    label 'with_gpus'

    publishDir path: "${params.outdir}/${params.sampleid}/results/megalodon_${task.process}/", mode: 'copy'

    input:
    path genomeref
    path reads
    // path variants

    output:
    path "megalodon_results"

    script:
    """
    export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    megalodon \
        --devices "cuda:all" \
        --processes 8 \
        --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
        --guppy-params "-d ${params.rerio_base}/basecall_models/" \
        --guppy-config ${params.megalodon_model} \
        --reference $genomeref \
        --outputs mod_mappings per_read_mods \
        --mappings-format bam \
        --mod-motif m CG 0 \
        --output-directory ./megalodon_results \
        $reads
        ## --outputs mod_mappings per_read_mods variant_mappings per_read_variants \
    """
}
        // --mod-motif Y A 0 \
        // --mod-motif Z CG 0 \



/* 
* Identification of modified bases using Megalodon
*/
process megalodon_aggregate {
    label 'megalodon'
    label 'bigmem'

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path megalodon_results
    path variants
    path genomeref

    output:
    path "mod_mappings_sort.bam", emit: mod_basecalls

    script:
    """
    megalodon_extras aggregate run \
    """
}

