
/* 
* Identification of modified bases using Megalodon
*/
process megalodon {
    label 'megalodon'
    label 'with_p100'

    publishDir path: "${params.outdir}/${params.sampleid}/results/megalodon_${task.process}/", mode: 'copy'

    input:
    path genomeref
    path reads
    // path variants

    output:
    path "megalodon_results", emit: megalodon_results

    script:
    """
    export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    megalodon \
        --devices "cuda:all" \
        --processes 18 \
        --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
        --guppy-params "-d ${params.rerio_base}/basecall_models/ --chunk_size 3000" \
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
process megalodon_withvariants {
    label 'megalodon'
    label 'gpu_debug'

    publishDir path: "${params.outdir}/${params.sampleid}/results/megalodon_${task.process}/", mode: 'copy'

    input:
    path genomeref
    path reads
    path variants

    output:
    path "megalodon_results"

    script:
    """
    export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    megalodon \
        --devices "cuda:all" \
        --processes 18 \
        --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
        --guppy-params "-d ${params.rerio_base}/basecall_models/" \
        --guppy-config ${params.megalodon_model} \
        --reference $genomeref \
        --variant-filename $variants \
        --outputs mod_mappings per_read_mods variant_mappings per_read_variants \
        --mappings-format bam \
        --mod-motif m CG 0 \
        --output-directory ./megalodon_results \
        $reads
    """
}


/* 
* Identification of modified bases using Megalodon
*/
process megalodon_aggregate {
    label 'megalodon'
    label 'cpu_high'
    label 'mem_mid'
    label 'time_high'

    publishDir path: "${params.outdir}/${params.sampleid}/results/${task.process}/", mode: 'copy'

    input:
    path megalodon_results
    path hap1ids
    path hap2ids

    output:
    path "megalodon_results/modified_bases.hap*.bed"

    script:
    """
    megalodon_extras aggregate run \
        --processes $task.cpus \
        --outputs mods \
        --megalodon-directory $megalodon_results \
        --read-ids-filename $hap1ids \
        --output-suffix hap1

    megalodon_extras aggregate run \
        --processes $task.cpus \
        --outputs mods \
        --megalodon-directory $megalodon_results \
        --read-ids-filename $hap2ids \
        --output-suffix hap2
    """
}

