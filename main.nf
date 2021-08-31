#!/usr/bin/env nextflow
/*
====================================
            N A N O W G S
====================================
 nanowgs Analysis Pipeline.
 https://github.com/jdemeul/nanowgs
------------------------------------
*/

nextflow.enable.dsl=2


// Command line shortcuts, quick entry point:

/* 
* Guppy basecalling
*/
workflow guppy_basecalling_cli {

    include { basecall_reads } from './modules/guppy'

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    basecall_reads( Channel.fromPath( params.ont_base_dir ), genomeref )

}


/* 
* Sam to sorted bam conversion and indexing
*/
workflow sam_to_sorted_bam {
    take:
        genomeref
        sam

    main:
        include { sam_to_sorted_bam as samtobam } from './modules/samtools'
        samtobam( sam, genomeref )

    emit:
        sorted_bam = samtobam.out.sorted_bam
        bam_index = samtobam.out.bam_index

}


/* 
* Sam to sorted bam conversion and indexing – CLI shortcut
*/
workflow sam_to_sorted_bam_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    sam = Channel.fromPath( params.mapped_sam, checkIfExists: true )

    sam_to_sorted_bam( sam, genomeref )

}


/* 
* Call structural variation and generate consensus
*/
workflow call_svs {
    take:
        genomeref
        bam
        bam_index

    main:
        include { sniffles_sv_calling as sniffles } from './modules/sniffles'
        include { svim_sv_calling as svim } from './modules/svim'
        include { cutesv_sv_calling as cutesv } from './modules/cutesv'
        include { svim_sv_filtering as filter } from './modules/bcftools'
        include { survivor_sv_consensus as survivor } from './modules/survivor'
    
        cutesv( bam, aligned_reads_idx, genomeref )
        sniffles( bam, bam_index )
        svim( bam, bam_index, genomeref )
        filter( svim.out.sv_calls )

        allsvs = cutesv.out.sv_calls
                    .mix( sniffles.out.sv_calls, svim.out.sv_calls_q10 )
                    .collect()
        survivor( allsvs )
    
    emit:
        cutesv = cutesv.out.sv_calls
        sniffles = sniffles.out.sv_calls
        svim = svim.out.sv_calls_q10
        consensus = survivor.out.sv_consensus

    
}


/* 
* Call structural variation and generate consensus
*/
workflow call_svs_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    aligned_reads = Channel.fromPath( params.aligned_bam )
    aligned_reads_idx = Channel.fromPath( params.aligned_bam + ".bai" )
    
    call_svs( genomeref, bam, bam_index )
    
}

/* 
* Call small variants using ONT Medaka
*/
workflow medaka_variant_calling_cli {

    include { medaka_snv_calling } from './modules/medaka'

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    aligned_reads = Channel.fromPath( params.aligned_bam )
    aligned_reads_idx = Channel.fromPath( params.aligned_bam + ".bai" )

    medaka_snv_calling( aligned_reads, aligned_reads_idx, genomeref )

}


/* 
* Call small variants using PEPPER-Margin-DeepVariant
*/
workflow pepper_deepvariant_calling_cli {

    include { deepvariant_snv_calling } from './modules/deepvariant'

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    aligned_reads = Channel.fromPath( params.aligned_bam )
    aligned_reads_idx = Channel.fromPath( params.aligned_bam + ".bai" )

    deepvariant_snv_calling( aligned_reads, aligned_reads_idx, genomeref )

}


/* 
* Run a de novo genome assembly using Shasta
*/
workflow shasta_assembly {
    take:
        fastq
        config
    
    main:
        include { run_shasta_assembly as shasta } from './modules/shasta'
        
        shasta( fastq, config )

    emit:
        shasta.out.assembly

}


/* 
* Run a de novo genome assembly using Shasta – CLI shortcut
*/
workflow shasta_assembly_cli {

    fastq = Channel.fromPath( params.processed_reads )
    config = Channel.fromPath( params.shasta_config )

    shasta_assembly( fastq, config )

}



/* 
* Process reads/squiggles from an ONT run
*/
workflow process_reads {
    take:
        genomeref
        ont_base

    main:
    
        // include { create_minimap_index } from './modules/minimap2'
        include { basecall_reads as basecall } from './modules/guppy'
        include { filter_reads as filter } from './modules/fastp'

        if ( params.rebasecall ) {

            // if ( !file( params.genomeref + "/indexes/minimap2-ont/genome.mmi" ).exists() ) {
            //     create_minimap_index( genomeref )
            //     genomeindex = create_minimap_index.out.mmi
            // } else {
            //     genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
            // }

            basecall( Channel.fromPath( params.ont_base_dir ), genomeref )
            filter( basecall.out.fastqs.collect() )
        } else {
            filter( Channel.fromPath( params.ont_base_dir + "**.fastq.gz" ).collect() )
        }

    emit:
        fastq_trimmed = filter.out.fastq_trimmed

}


/* 
* Process reads/squiggles from an ONT run – CLI shortcut
*/
workflow process_reads_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    ont_base = Channel.fromPath( params.ont_base_dir )

    process_reads( genomeref, ont_base )

}


/* 
* Align reads to a reference genome using minimap2 and turn into sorted bam
*/
workflow minimap_alignment {
    take: 
        fastqs
        genomeref

    main:
        include { minimap_alignment as minimap } from './modules/minimap2'
        // include { sam_to_sorted_bam as samtobam } from './modules/samtools'

        genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
        fastqs = Channel.fromPath( params.processed_reads )
    
        // use index if matched index is available, otherwise do on the fly
        genomeindex = file( file(params.genomeref).getParent() + "/" + file(params.genomeref).getSimpleName() + ".mmi" )
        if ( genomeindex.exists() ) {
            minimap( Channel.fromPath( genomeindex ), fastqs )
        } else {
            minimap( genomeref, fastqs )
        }

        // alignment and conversion into indexed sorted bam
        // samtobam( minimap.out.mapped_sam, genomeref )

    emit:
        sam = minimap.out.mapped_sam
        // sorted_bam = samtobam.out.sorted_bam
        // bam_index = samtobam.out.bam_index

}


/* 
* Align reads to a reference genome using minimap2 and turn into sorted bam – CLI shortcut
*/
workflow minimap_alignment_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    fastqs = Channel.fromPath( params.processed_reads )

    minimap_alignment( fastqs, genomeref )

}



workflow lra_alignment_sv_calling {
    take: 
        fastqs
        genomeref

    main:
        include { create_lra_index; lra_alignment } from './modules/lra'
        include { sam_to_sorted_bam } from './modules/samtools'
        include { cutesv_sv_calling } from './modules/cutesv'

        // genome indexing
        if ( !file( params.genomeref + "/indexes/lra-ont/genome.fa.gli" ).exists() || !file( params.genomeref + "/indexes/lra-ont/genome.fa.mmi" ).exists() ) {
            create_lra_index( genomeref )
            genomeindex_gli = create_lra_index.out.gli
            genomeindex_mmi = create_lra_index.out.mmi
        } else {
            genomeindex_gli = Channel.fromPath( params.genomeref + "/indexes/lra-ont/genome.fa.gli" )
            genomeindex_mmi = Channel.fromPath( params.genomeref + "/indexes/lra-ont/genome.fa.mmi" )
        }

        // alignment and conversion into indexed sorted bam
        lra_alignment( genomeref, fastqs, genomeindex_gli, genomeindex_mmi )
        sam_to_sorted_bam( lra_alignment.out.mapped_sam, genomeref, lra_alignment.out.aligner )

        // SV calling using cuteSV
        cutesv_sv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

}






workflow {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    ont_base = Channel.fromPath( params.ont_base_dir, checkIfExists: true )
    // reads = Channel.fromPath( params.processed_reads, checkIfExists: true )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    // process reads
    process_reads( genomeref, ont_base )

    // 1. start assembly
    shasta_assembly( process_reads.out.fastq_trimmed )

    // 2. reference alignment and conversion into indexed sorted bam
    minimap_alignment( genomeref, process_reads.out.fastq_trimmed )
    sam_to_sorted_bam( minimap_alignment.out.sam, genomeref )

    sv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

    // lra_alignment_sv_calling( process_reads.out.fastq_trimmed, genomeref )
    // minimap_alignment_snv_calling( process_reads.out.fastq_trimmed, genomeref )

}

// -m ${Math.floor( $task.memory / $task.cpus ) }
// ch_reference_fasta.view()
// ch_reference_index.view()