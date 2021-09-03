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
include { survivor_sv_consensus as survivor } from './modules/survivor'

include { medaka_snv_calling as medaka_snv } from './modules/medaka'
include { deepvariant_snv_calling as deepvariant } from './modules/deepvariant'

include { run_shasta_assembly as shasta } from './modules/shasta'
include { racon_assembly_polishing as racon } from './modules/racon'
include { medaka_assembly_polishing as medaka_polish } from './modules/medaka'
include { medaka_assembly_polish_align; medaka_assembly_polish_stitch; medaka_assembly_polish_consensus } from './modules/medaka'

include { create_lra_index; lra_alignment } from './modules/lra'



/* 
* Guppy basecalling
*/
workflow guppy_basecalling_cli {

    // include { basecall_reads } from './modules/guppy'

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    basecall( Channel.fromPath( params.ont_base_dir ), genomeref )

}


// /* 
// * Sam to sorted bam conversion and indexing
// */
// workflow sam_to_sorted_bam {
//     take:
//         genomeref
//         sam

//     main:
//         // include { sam_to_sorted_bam as samtobam } from './modules/samtools'
//         samtobam( sam, genomeref )

//     emit:
//         sorted_bam = samtobam.out.sorted_bam
//         bam_index = samtobam.out.bam_index

// }


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

    main:
        // include { sniffles_sv_calling as sniffles } from './modules/sniffles'
        // include { svim_sv_calling as svim } from './modules/svim'
        // include { cutesv_sv_calling as cutesv } from './modules/cutesv'
        // include { svim_sv_filtering as filtersvim } from './modules/bcftools'
        // include { survivor_sv_consensus as survivor } from './modules/survivor'
    
        cutesv( bam, bam_index, genomeref )
        sniffles( bam, bam_index )
        svim( bam, bam_index, genomeref )
        filtersvim( svim.out.sv_calls )

        allsvs = cutesv.out.sv_calls
                    .mix( sniffles.out.sv_calls, filtersvim.out.sv_calls_q10 )
                    .collect()
        survivor( allsvs )
    
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
    
    call_svs( genomeref, bam, bam_index )
    
}

/* 
* Call small variants using ONT Medaka
*/
workflow medaka_variant_calling_cli {

    // include { medaka_snv_calling } from './modules/medaka'

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    aligned_reads = Channel.fromPath( params.aligned_bam )
    aligned_reads_idx = Channel.fromPath( params.aligned_bam + ".bai" )

    medaka_snv( aligned_reads, aligned_reads_idx, genomeref )

}


// /* 
// * Call small variants using PEPPER-Margin-DeepVariant
// */
// workflow pepper_deepvariant_calling {
//     take:
//         bam
//         bam_index
//         genomeref
    
//     main:
//         // include { deepvariant_snv_calling as deepvariant } from './modules/deepvariant'
//         deepvariant( bam, bam_index, genomeref )

//     emit:
//         snvs = deepvariant.out.indel_snv_vcf
// }


/* 
* Call small variants using PEPPER-Margin-DeepVariant – CLI shortcut
*/
workflow pepper_deepvariant_calling_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    bam = Channel.fromPath( params.aligned_bam )
    bam_index = Channel.fromPath( params.aligned_bam + ".bai" )

    deepvariant( bam, bam_index, genomeref )

}


// /* 
// * Run a de novo genome assembly using Shasta
// */
// workflow shasta_assembly {
//     take:
//         fastq
//         config
    
//     main:
//         // include { run_shasta_assembly as shasta } from './modules/shasta'
        
//         shasta( fastq, config )

//     emit:
//         assembly = shasta.out.assembly

// }


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
        // include { racon_assembly_polishing as racon } from './modules/racon'
        // include { medaka_assembly_polishing as medaka } from './modules/medaka'
        // include { minimap_alignment as minimap } from './modules/minimap2'

        minimap( assembly, fastq )
        racon( fastq, minimap.out.mapped_sam, assembly )
        // medaka_polish( fastq, racon.out.consensus )
        medaka_polish_parallel( fastq, racon.out.consensus )

    emit:
        polished_assembly = medaka_polish_parallel.out.consensus
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
        // include { basecall_reads as basecall } from './modules/guppy'
        // include { filter_reads as filter } from './modules/fastp'

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
// workflow minimap_alignment {
//     take: 
//         genomeref
//         fastqs

//     main:
//         // include { minimap_alignment as minimap } from './modules/minimap2'
//         // include { sam_to_sorted_bam as samtobam } from './modules/samtools'
    
//         // use index if matched index is available, otherwise do on the fly
//         genomeindex = file( file(params.genomeref).getParent() + "/" + file(params.genomeref).getSimpleName() + ".mmi" )
//         if ( genomeindex.exists() ) {
//             minimap( Channel.fromPath( genomeindex ), fastqs )
//         } else {
//             minimap( genomeref, fastqs )
//         }

//         // alignment and conversion into indexed sorted bam
//         // samtobam( minimap.out.mapped_sam, genomeref )

//     emit:
//         sam = minimap.out.mapped_sam
//         // sorted_bam = samtobam.out.sorted_bam
//         // bam_index = samtobam.out.bam_index

// }


/* 
* Align reads to a reference genome using minimap2 and turn into sorted bam – CLI shortcut
*/
workflow minimap_alignment_cli {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )
    fastqs = Channel.fromPath( params.processed_reads )

    minimap( genomeref, fastqs )

}



workflow lra_alignment_sv_calling {
    take: 
        fastqs
        genomeref

    main:
        // include { create_lra_index; lra_alignment } from './modules/lra'
        // include { sam_to_sorted_bam } from './modules/samtools'
        // include { cutesv_sv_calling } from './modules/cutesv'

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
        assembly_polishing( shasta.out.assembly, fastq)

        minimap( assembly_polishing.out.polished_assembly, fastq )
        samtobam( minimap.out.mapped_sam, assembly_polishing.out.polished_assembly )

        call_svs( assembly_polishing.out.polished_assembly, samtobam.out.sorted_bam, samtobam.out.bam_index )
        deepvariant( samtobam.out.sorted_bam, samtobam.out.bam_index, assembly_polishing.out.polished_assembly )

    emit:
        polished_assembly = assembly_polishing.out.polished_assembly
        svs = call_svs.out.consensus
        snvs = deepvariant.out.indel_snv_vcf

}


workflow reference_based_variant_calling {
    take:
        fastq
        genomeref
    
    main:
        minimap( genomeref, fastq )
        samtobam( minimap.out.mapped_sam, genomeref )

        call_svs( genomeref, samtobam.out.sorted_bam, samtobam.out.bam_index )
        deepvariant( samtobam.out.sorted_bam, samtobam.out.bam_index, genomeref )

    emit:
        svs = call_svs.out.consensus
        snvs = deepvariant.out.indel_snv_vcf

}


workflow medaka_polish_parallel {
    take:
        fastqs
        draft

    main:
        medaka_assembly_polish_align( fastqs, draft )
        // draft.splitFasta( record: [id: true, seqString: false ]).view { it.id }
        contigs = draft.splitFasta( record: [id: true, seqString: false ]).map { it.id }

        medaka_assembly_polish_consensus( medaka_assembly_polish_align.out.calls_to_draft,
                                         medaka_assembly_polish_align.out.calls_to_draft_index,
                                         contigs )

        medaka_assembly_polish_stitch( medaka_assembly_polish_consensus.out.probs.collect() )
    
    emit:
        conensus = medaka_assembly_polish_stitch.out.consensus
}

workflow medaka_polish_parallel_cli {

    fastqs = Channel.fromPath( params.processed_reads )
    draft = Channel.fromPath( params.genomeref )

    medaka_assembly_polish_parallel( fastqs, draft )    
}


// workflow trial_cli {

//         // medaka_assembly_polishing_align( fastqs, draft )
//         draft = Channel.fromPath( params.genomeref )
//         draft.splitFasta( record: [id: true, seqString: false ]).map { it.id }

//         // medaka_assembly_polishing_align( medaka_assembly_polishing_consensus.out.bam,
//         //                                  medaka_assembly_polishing_consensus.out.bamindex,
//         //                                  draft.splitFasta( record: [id: true, seqString: false ]) )
    
// }

workflow {

    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    ont_base = Channel.fromPath( params.ont_base_dir, checkIfExists: true )
    // reads = Channel.fromPath( params.processed_reads, checkIfExists: true )
    // genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    // process reads
    process_reads( genomeref, ont_base )

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