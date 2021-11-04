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

include { basecall_reads as basecall } from './modules/guppy'
include { filter_reads as filter } from './modules/fastp'

include { minimap_alignment as minimap } from './modules/minimap2'
include { sam_to_sorted_bam as samtobam } from './modules/samtools'

include { sniffles_sv_calling as sniffles } from './modules/sniffles'
include { svim_sv_calling as svim } from './modules/svim'
include { cutesv_sv_calling as cutesv } from './modules/cutesv'
include { svim_sv_filtering as filtersvim } from './modules/bcftools'
include { vcf_concat } from './modules/bcftools'

include { survivor_sv_consensus as survivor } from './modules/survivor'

include { megalodon } from './modules/megalodon'

include { medaka_snv_calling as medaka_snv } from './modules/medaka'
include { deepvariant_snv_calling as deepvariant } from './modules/deepvariant'
// include { deepvariant_snv_calling_gpu_parallel as deepvariant_par } from './modules/deepvariant'

include { run_shasta_assembly as shasta } from './modules/shasta'
include { racon_assembly_polishing as racon } from './modules/racon'
include { medaka_assembly_polishing as medaka_polish } from './modules/medaka'
include { medaka_assembly_polish_align; medaka_assembly_polish_stitch; medaka_assembly_polish_consensus } from './modules/medaka'

include { create_lra_index; lra_alignment } from './modules/lra'



/* 
* Guppy basecalling
*/
workflow guppy_basecalling_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    basecall( Channel.fromPath( params.ont_base_dir ), genomeref )

}


/* 
* Sam to sorted bam conversion and indexing – CLI shortcut
*/
workflow sam_to_sorted_bam_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    sam = Channel.fromPath( params.mapped_sam, checkIfExists: true )

    samtobam( sam, genomeref )

}


/* 
* Call structural variation and generate consensus
*/
workflow call_svs {
    take:
        genomeref
        bam
        bam_index
        step

    main:    
        cutesv( bam, bam_index, genomeref, step )
        sniffles( bam, bam_index, step )
        svim( bam, bam_index, genomeref, step )
        filtersvim( svim.out.sv_calls, step )

        allsvs = cutesv.out.sv_calls
                    .mix( sniffles.out.sv_calls, filtersvim.out.sv_calls_q10 )
                    .collect()
        survivor( allsvs, step )
    
    emit:
        cutesv = cutesv.out.sv_calls
        sniffles = sniffles.out.sv_calls
        svim = filtersvim.out.sv_calls_q10
        consensus = survivor.out.sv_consensus

    
}


/* 
* Call structural variation and generate consensus
*/
workflow call_svs_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    bam = Channel.fromPath( params.aligned_bam )
    bam_index = Channel.fromPath( params.aligned_bam + ".bai" )
    
    call_svs( genomeref, bam, bam_index, "cli" )
    
}

/* 
* Call small variants using ONT Medaka
*/
workflow medaka_variant_calling_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    aligned_reads = Channel.fromPath( params.aligned_bam )
    aligned_reads_idx = Channel.fromPath( params.aligned_bam + ".bai" )

    medaka_snv( aligned_reads, aligned_reads_idx, genomeref )

}


/* 
* Call small variants using PEPPER-Margin-DeepVariant – CLI shortcut
*/
workflow pepper_deepvariant_calling_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    bam = Channel.fromPath( params.aligned_bam )
    bam_index = Channel.fromPath( params.aligned_bam + ".bai" )

    deepvariant( bam, bam_index, genomeref, "cli" )

}


/* 
* Run a de novo genome assembly using Shasta – CLI shortcut
*/
workflow shasta_assembly_cli {

    fastq = Channel.fromPath( params.processed_reads )
    config = Channel.fromPath( params.shasta_config )

    shasta( fastq, config )

}


workflow assembly_polishing {
    take:
        assembly
        fastq
    
    main:
        minimap( assembly, fastq )
        racon( fastq, minimap.out.mapped_sam, assembly )
        // medaka_polish( fastq, racon.out.consensus )
        medaka_polish_parallel( fastq, racon.out.consensus )

    emit:
        polished_assembly = medaka_polish_parallel.out.consensus
}


workflow deepvariant_parallel {
    take:
        bam
        bam_index
        genomeref

    main:
        contigs = genomeref.splitFasta( record: [id: true, seqString: false ]).map { it.id }

        deepvariant_par( bam, bam_index, genomeref, contigs )

        vcf_concat( deepvariant_par.out.indel_snv_vcf.collect() )
    
    emit:
        indel_snv_vcf = vcf_concat.out.merged_vcf
}


workflow medaka_polish_parallel {
    take:
        fastqs
        draft

    main:
        medaka_assembly_polish_align( fastqs, draft )
        // draft.splitFasta( record: [id: true, seqString: false ]).view { it.id }
        contigs = draft.splitFasta( record: [id: true, seqString: false ]).map { it.id }.randomSample( 100000, 234 ).buffer( size: 50, remainder: true )

        medaka_assembly_polish_consensus( medaka_assembly_polish_align.out.calls_to_draft,
                                         medaka_assembly_polish_align.out.calls_to_draft_index,
                                         contigs )

        medaka_assembly_polish_stitch( medaka_assembly_polish_consensus.out.probs.collect(), draft )
    
    emit:
        consensus = medaka_assembly_polish_stitch.out.consensus
}

workflow medaka_polish_parallel_cli {

    fastqs = Channel.fromPath( params.processed_reads )
    draft = Channel.fromPath( params.genomeref )

    medaka_assembly_polish_parallel( fastqs, draft )    
}



/* 
* Process reads/squiggles from an ONT run
*/
workflow process_reads {
    take:
        genomeref
        ont_base

    main:

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
* Align reads to a reference genome using minimap2 and turn into sorted bam – CLI shortcut
*/
workflow minimap_alignment_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    fastqs = Channel.fromPath( params.processed_reads )

    minimap( genomeref, fastqs )

}


workflow minimap_align_bamout {
    take:
        genomeref
        fastq
    
    main:
        minimap( genomeref, fastq )
        samtobam( minimap.out.mapped_sam, genomeref )

    emit:
        bam = samtobam.out.sorted_bam
        idx = samtobam.out.bam_index

}


workflow lra_alignment_sv_calling {
    take: 
        fastqs
        genomeref

    main:
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
        samtobam( lra_alignment.out.mapped_sam, genomeref )

        // SV calling using cuteSV
        cutesv( samtobam.out.sorted_bam, samtobam.out.bam_index, genomeref )

}


workflow assembly_based_variant_calling {
    take:
        fastq
        // genomeref
    
    main:
        // 1. Assembly-based pipeline
        shasta( fastq, params.shasta_config )
        // assembly_polishing( shasta.out.assembly, fastq)

        // minimap_align_bamout( assembly_polishing.out.polished_assembly, fastq )
    //     minimap( assembly_polishing.out.polished_assembly, fastq )
    //     samtobam( minimap.out.mapped_sam, assembly_polishing.out.polished_assembly )

        // call_svs( assembly_polishing.out.polished_assembly, minimap_align_bamout.out.bam, minimap_align_bamout.out.idx, "assembly" )
        // deepvariant( minimap_align_bamout.out.bam, minimap_align_bamout.out.idx, assembly_polishing.out.polished_assembly, "assembly" )
    //     call_svs( assembly_polishing.out.polished_assembly, samtobam.out.sorted_bam, samtobam.out.bam_index )
    //     deepvariant( samtobam.out.sorted_bam, samtobam.out.bam_index, assembly_polishing.out.polished_assembly )

    emit:
        polished_assembly = shasta.out.assembly
        // svs = call_svs.out.consensus
        // snvs = deepvariant.out.indel_snv_vcf
        // snvs_idx = deepvariant.out.indel_snv_vcf_index

}


workflow reference_based_variant_calling {
    take:
        fastq
        genomeref
    
    main:
        // minimap( genomeref, fastq )
        // samtobam( minimap.out.mapped_sam, genomeref )
        minimap_align_bamout( genomeref, fastq )

        call_svs( genomeref, minimap_align_bamout.out.bam, minimap_align_bamout.out.idx, "reference" )

        // note that on the current VSC system (P100 and V100 nodes) GPU nodes
        // lack sufficient memory to run PEPPER-DeepVariant genome-wide
        // the current if statement reflects that and splits up the genome by chromosome to run as separate jobs
        // if ( params.deepvariant_with_gpu ) {
            // deepvar = deepvariant_parallel( minimap_align_bamout.out.bam, minimap_align_bamout.out.idx, genomeref )
        // } else {
            deepvariant( minimap_align_bamout.out.bam, minimap_align_bamout.out.idx, genomeref, "reference" )
        // }

    emit:
        svs = call_svs.out.consensus
        snvs = deepvariant.out.indel_snv_vcf
        snvs_idx = deepvariant.out.indel_snv_vcf_index

}


workflow {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    ont_base = Channel.fromPath( params.ont_base_dir, checkIfExists: true )
    // reads = Channel.fromPath( params.processed_reads, checkIfExists: true )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    // process reads
    process_reads( genomeref, ont_base )

    // start megalodon
    // megalodon( genomeref, ont_base )

    // assembly based variant calling
    assembly_based_variant_calling( process_reads.out.fastq_trimmed )

    // 2. Reference alignment-based pipeline
    reference_based_variant_calling( process_reads.out.fastq_trimmed, genomeref )

    // lra_alignment_sv_calling( process_reads.out.fastq_trimmed, genomeref )
    // minimap_alignment_snv_calling( process_reads.out.fastq_trimmed, genomeref )

}

// -m ${Math.floor( $task.memory / $task.cpus ) }
// ch_reference_fasta.view()
// ch_reference_index.view()