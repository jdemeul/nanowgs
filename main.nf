#!/usr/bin/env nextflow
/*
====================================
            N A N O W G S
====================================
 nanowgs Analysis Pipeline.
 https://github.com/jdemeul/nanowgs
------------------------------------
*/


/* 
 * pipeline input parameters 
 */
// params {
params.genomeref           = "/staging/leuven/stg_00002/lcb/jdemeul/reference/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta"
params.guppy_gpu           = true
params.outdir              = "./results"
params.tracedir            = "${params.outdir}/pipeline_info"
// }


/* 
 * index a reference genome with minimap2
 */
process minimap_index {
    // tag "$genomeref"
    // label 'process_medium'

    input:
    path genomeref from params.genomeref

    output:
    path 'index' into ch_index

    publishDir path: '.', mode: 'copy',
        saveAs: { "${file(params.genomeref).getParent()}/${file(params.genomeref).getBaseName()}_map-ont.mmi"}

    script:
    """
    minimap2 -ax map-ont -t $task.cpus -d index $genomeref
    """
}

ch_index.view()