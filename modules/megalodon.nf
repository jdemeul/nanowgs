
/* 
* Identification of modified bases using Megalodon
*/
process megalodon_modifications {
    label 'megalodon'
    label ( params.with_gpu ? 'with_gpu': 'cpu_high')
    label ( params.with_gpu ? null: 'mem_high')
    label ( params.with_gpu ? null: 'time_high')

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path reads
    path variants
    path genomeref

    output:
    path "mod_mappings_sort.bam", emit: mod_basecalls

    script:
    if ( ${params.with_gpu} ) 
        """
        export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
        megalodon \
            --devices "cuda:all" \
            --processes $task.cpus \
            --guppy-server-path /opt/ont/ont-guppy/bin/guppy_basecall_server \
            --guppy-params "-d /staging/leuven/stg_00002/lcb/jdemeul/software/rerio/basecall_models/" \
            --guppy-config ${params.megalodon_model} \
            --reference $genomeref \
            --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mods \
            --mappings-format bam \
            --overwrite \
            --sort-mappings \
            --write-mods-text \
            --mod-motif Y A 0 \
            --mod-motif Z CG 0 \
            --output-directory /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/results/20210328_MouseBrain_Hia5_megalodon/ \
            /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/data/20210326_1208_MN34250_FAP14858_0224c788/
        """
    else 
        """
        megalodon \
            --processes $task.cpus \
            --guppy-server-path /opt/ont/ont-guppy/bin/guppy_basecall_server \
            --guppy-params "-d /staging/leuven/stg_00002/lcb/jdemeul/software/rerio/basecall_models/" \
            --guppy-config ${params.megalodon_model} \
            --reference $genomeref \
            --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mods \
            --mappings-format bam \
            --overwrite \
            --sort-mappings \
            --write-mods-text \
            --mod-motif Y A 0 \
            --mod-motif Z CG 0 \
            --output-directory /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/results/20210328_MouseBrain_Hia5_megalodon/ \
            /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/data/20210326_1208_MN34250_FAP14858_0224c788/
        """


}

