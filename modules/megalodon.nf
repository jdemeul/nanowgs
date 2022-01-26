
/* 
* Identification of modified bases using Megalodon
*/
process megalodon {
    label 'megalodon'
    label 'with_p100node'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path genomeref
    path reads
    // path variants

    output:
    path "megalodon_results", emit: megalodon_results

    script:

    if ( params.megalodon_modmotif2 != "none" )
        """
        # export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
        megalodon \
            --devices ${params.gpu_devices} \
            --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
            --guppy-params "-d /rerio/basecall_models/" \
            --guppy-config ${params.megalodon_model} \
            --reference $genomeref \
            --outputs basecalls mod_mappings per_read_mods \
            --mappings-format bam \
            --mod-motif ${params.megalodon_modmotif} \
            --mod-motif ${params.megalodon_modmotif2} \
            --processes 16 \
            --guppy-concurrent-reads 40 \
            --guppy-timeout 120 \
            --output-directory ./megalodon_results \
            --num-read-enumeration-threads 1 \
            --num-extract-signal-processes 2 \
            $reads
            ## --outputs mod_mappings per_read_mods variant_mappings per_read_variants --chunk_size 3000 \
        """
    else
        """
        # export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
        megalodon \
            --devices ${params.gpu_devices} \
            --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server \
            --guppy-params "-d /rerio/basecall_models/" \
            --guppy-config ${params.megalodon_model} \
            --reference $genomeref \
            --outputs basecalls mod_mappings per_read_mods \
            --mappings-format bam \
            --mod-motif ${params.megalodon_modmotif} \
            --processes 16 \
            --guppy-concurrent-reads 40 \
            --guppy-timeout 120 \
            --output-directory ./megalodon_results \
            --num-read-enumeration-threads 1 \
            --num-extract-signal-processes 2 \
            $reads
        """
}
        // --mod-motif Y A 0 \
        // --mod-motif Z CG 0 \
// node optimisation on P100
// singularity run --nv -B /staging/leuven/stg_00002/lcb/gc_test/ -B /staging/leuven/stg_00002/lcb/jdemeul/ /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-guppy-5.0.16-megalodon_v2.3.5-6d32ce0-rerio.img megalodon         --devices "cuda:all"         --guppy-server-path /opt/ont-guppy/bin/guppy_basecall_server         --guppy-params "-d /rerio/basecall_models/"         --guppy-config res_dna_r941_min_modbases_5mC_v001.cfg         --reference /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/indexes/minimap2-ont/genome.fa --outputs basecalls mod_mappings per_read_mods         --mappings-format bam         --mod-motif m CG 0         --processes 16         --guppy-concurrent-reads 40         --guppy-timeout 120    --overwrite     --output-directory ./megalodon_results --num-read-enumeration-threads 1 --num-extract-signal-processes 1 /staging/leuven/stg_00002/lcb/gc_test/20210719_ASA_Edin_BA24_14_18/20210719_1840_3A_PAH21935_2e6023c0/fast5_pass/

// # 2x P100 GPU
// 16.4/1.94e6 for 4 proc and 2000 chunks
// 11.62/1.4e6 for 4 proc 10 read 3000 chunks
// 30.6/3.59e6 for 4 proc 20 reads
// 30.89/3.59e6 for 8 proc 10 reads
// 62/7.5e6 for 8 proc 20 reads
// 82/1e7 for 16 proc 20 reads
// 72/9.1e6 for 10 proc 20 reads
// 81/1e7 for 12 proc 25 reads
// 81/1e7 for 10 proc 30 reads
// 79/1e7 for 6 proc 50 reads
// 81/1e7 for 8 proc 40 reads

// # 4x P100 GPU
// 125/1.6e7 for 8 proc 40 reads
// 160/2e7 for 16 proc 40 reads



/* 
* Identification of modified bases using Megalodon
*/
process megalodon_withvariants {
    label 'megalodon'
    label 'gpu_debug'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

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
    label 'mem_high'
    label 'time_high'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

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

