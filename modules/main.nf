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



// /* 
// * SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
// */
// process create_personal_genome {
//     label 'process_low'
//     label 'deepvariant'

//     publishDir path: "${params.outdir}/results/crossstitch", mode: 'copy'

//     input:
//     path phased_snps
//     path unphased_svs
//     path sorted_bam
//     path bam_index
//     path genomeref
//     val karyotype
//     val refine

//     output:
//     path "*.vcf.gz", emit: indel_snv_vcf
//     path "*.vcf.gz.tbi", emit: indel_snv_vcf_index
//     path "logs"
//     path "intermediate_files"
//     path "*.visual_report.html"

//     script:
//     """
//     export XDG_CONFIG_HOME="/staging/leuven/stg_00002/lcb/jdemeul/software/"
//     module use /staging/leuven/stg_00002/lcb/jdemeul/software/easybuild/modules/all/
//     ml --ignore-cache Java/13.0.2 SAMtools/1.9-GCC-6.4.0-2.28 HTSlib/1.9-GCC-6.4.0-2.28 pigz/2.6-GCCcore-6.4.0 Racon/1.4.13-GCCcore-9.3.0
//     /staging/leuven/stg_00002/lcb/jdemeul/software/crossstitch/src/crossstitch.sh  \
//         /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20.phased_PASS.vcf /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20_LRA_cuteSV_svs_PASS.vcf /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20_minimap2.bam /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chr20.fasta ASA_Edin_BA24_14_18_chr20 xy 1
//     run_pepper_margin_deepvariant call_variant \
//         -b $sorted_bam \
//         -f $genomeref \
//         -o . \
//         -p ${params.sampleid} \
//         -t $task.cpus \
//         --ont \
//         --phased_output
//     """
// }


// /* 
// * SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
// */
// process megalodon_modifications {
//     label 'process_high'
//     label 'megalodon'
//     label (${params.with_gpu} ? 'with_gpu': null)

//     publishDir path: "${params.outdir}/results/", mode: 'copy'

//     input:
//     path reads
//     path variants
//     path genomeref

//     output:
//     path "mod_mappings_sort.bam", emit: mod_basecalls

//     script:
//     if ( ${params.with_gpu} ) 
//         """
//         export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
//         megalodon \
//             --devices "cuda:all" \
//             --processes $task.cpus \
//             --guppy-server-path /opt/ont/ont-guppy/bin/guppy_basecall_server \
//             --guppy-params "-d /staging/leuven/stg_00002/lcb/jdemeul/software/rerio/basecall_models/" \
//             --guppy-config ${params.megalodon_model} \
//             --reference $genomeref \
//             --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mods \
//             --mappings-format bam \
//             --overwrite \
//             --sort-mappings \
//             --write-mods-text \
//             --mod-motif Y A 0 \
//             --mod-motif Z CG 0 \
//             --output-directory /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/results/20210328_MouseBrain_Hia5_megalodon/ \
//             /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/data/20210326_1208_MN34250_FAP14858_0224c788/
//         """
//     else 
//         """
//         megalodon \
//             --processes $task.cpus \
//             --guppy-server-path /opt/ont/ont-guppy/bin/guppy_basecall_server \
//             --guppy-params "-d /staging/leuven/stg_00002/lcb/jdemeul/software/rerio/basecall_models/" \
//             --guppy-config ${params.megalodon_model} \
//             --reference $genomeref \
//             --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mods \
//             --mappings-format bam \
//             --overwrite \
//             --sort-mappings \
//             --write-mods-text \
//             --mod-motif Y A 0 \
//             --mod-motif Z CG 0 \
//             --output-directory /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/results/20210328_MouseBrain_Hia5_megalodon/ \
//             /staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/data/20210326_1208_MN34250_FAP14858_0224c788/
//         """


// }


// /* 
// * Run de novo genome assembly using Shasta
// */
// process run_shasta_assembly {
//     label 'bigmem'
//     label 'shasta'

//     input:
//     path phased_snps

//     output:
//     path "*.vcf.gz", emit: indel_snv_vcf

//     script:
//     """
//     singularity run \
//     -B /staging/leuven/stg_00002/lcb/jdemeul/ \
//     -B /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/shasta_assembly/:/output \
//     /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-shasta-docker-latest.img 0.7.0 \
//     --input /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ASA_Edin_BA24_38_17_trimmed.fastq \
//     --conf /staging/leuven/stg_00002/lcb/jdemeul/software/shasta/conf/Nanopore-Sep2020.conf

//     singularity run -B /staging/leuven/stg_00002/lcb/jdemeul/ /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-mummer-4.0.0rc1.img nucmer -t 16 --maxmatch -l 100 -c 1000 -p /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/fasta/chm13_v1.1_chrY_KI270740_EBV.fasta /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/Assembly.fasta
//     python DotPrep.py --delta /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ASA_Edin_BA24_38_17_denovo_vs_CHM13v1.1.delta --out /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ASA_Edin_BA24_38_17_denovo_vs_CHM13v1.1.DotPrep.out
    
//     singularity run -B /staging/leuven/stg_00002/lcb/jdemeul/ /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-quast-5.1.0rc1.img quast.py \
//     /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/Assembly.fasta \
//     -r /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/fasta/chm13_v1.1_chrY_KI270740_EBV.fasta \
//     -g gene:/staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/annotation/chm13.draft_v1.1.gene_annotation.v4.gff3.gz \
//     -t 16 \
//     --large \
//     --k-mer-stats \
//     --ref-bam /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/bam/ASA_Edin_BA24_38_17_minimap2.bam \
//     -o /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/quast_output
//     """
// }


// }

workflow minimap_alignment_snv_calling {
    take: 
        fastqs
        genomeref
    main:
        // genomeref = Channel.fromPath( params.genomeref + "/fasta/genome.fa", checkIfExists: true  )
        // fastqs = Channel.fromPath( params.fastqs )
        
        // genome indexing
        if ( !file( params.genomeref + "/indexes/minimap2-ont/genome.mmi" ).exists() ) {
            create_minimap_index( genomeref )
            genomeindex = create_minimap_index.out.mmi
        } else {
            genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
        }

        // alignment and covnersion into indexed sorted bam
        minimap_alignment( genomeref, genomeindex, fastqs )
        sam_to_sorted_bam( minimap_alignment.out.mapped_sam, genomeref, minimap_alignment.out.aligner )

        if ( params.minimap_sv_calls ) {
            cutesv_sv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )
        }

        // SNV calling using PEPPER-margin-DeepVariant
        deepvariant_snv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

}


workflow sv_calling {
    take:
        bam
        bam_index
        genomeref
    main:

        cutesv_sv_calling( bam, bam_index, genomeref )
        sniffles_sv_calling( bam, bam_index )
        svim_sv_calling( bam, bam_index, genomeref )
        svim_sv_filtering( svim_sv_calling.out.sv_calls )
    
    emit:
        sv_cutesv = cutesv_sv_calling.out.sv_calls
        sv_svim = svim_sv_filtering.out.sv_calls_q10
        sv_sniffles = sniffles_sv_calling.out.sv_calls
}


// workflow lra_alignment_sv_calling {
//     take: 
//         fastqs
//         genomeref

//     main:
//         include { create_lra_index, lra_alignment } from './modules/lra'
//         include { sam_to_sorted_bam } from './modules/samtools'
//         include { cutesv_sv_calling } from './modules/cutesv'

//         // genome indexing
//         if ( !file( params.genomeref + "/indexes/lra-ont/genome.fa.gli" ).exists() || !file( params.genomeref + "/indexes/lra-ont/genome.fa.mmi" ).exists() ) {
//             create_lra_index( genomeref )
//             genomeindex_gli = create_lra_index.out.gli
//             genomeindex_mmi = create_lra_index.out.mmi
//         } else {
//             genomeindex_gli = Channel.fromPath( params.genomeref + "/indexes/lra-ont/genome.fa.gli" )
//             genomeindex_mmi = Channel.fromPath( params.genomeref + "/indexes/lra-ont/genome.fa.mmi" )
//         }

//         // alignment and conversion into indexed sorted bam
//         lra_alignment( genomeref, fastqs, genomeindex_gli, genomeindex_mmi )
//         sam_to_sorted_bam( lra_alignment.out.mapped_sam, genomeref, lra_alignment.out.aligner )

//         // SV calling using cuteSV
//         cutesv_sv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

// }


workflow process_reads {
    take:
        genomeref
        ont_base

    main:
    
        include { create_minimap_index } from './modules/minimap2'
        include { basecall_reads } from './modules/guppy'
        include { filter_reads } from './modules/fastp'

        if ( params.rebasecall ) {

            if ( !file( params.genomeref + "/indexes/minimap2-ont/genome.mmi" ).exists() ) {
                create_minimap_index( genomeref )
                genomeindex = create_minimap_index.out.mmi
            } else {
                genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
            }

            basecall_reads( Channel.fromPath( params.ont_base_dir ), genomeref, genomeindex )
            filter_reads( basecall_reads.out.fastqs.collect() )
        } else {
            filter_reads( Channel.fromPath( params.ont_base_dir + "**.fastq.gz" ).collect() )
        }

    emit:
        fastq_trimmed = filter_reads.out.fastq_trimmed

}


workflow guppy_basecalling {

    include { basecall_reads } from './modules/guppy'

    genomeref = Channel.fromPath( params.genomeref + "/fasta/genome.fa", checkIfExists: true  )
    genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )

    basecall_reads( Channel.fromPath( params.ont_base_dir ), genomeref, genomeindex )

}


workflow medaka_variant_calling {

    include { medaka_snv_calling } from './modules/medaka'

    genomeref = Channel.fromPath( params.genomeref + "/fasta/genome.fa", checkIfExists: true  )
    genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    aligned_reads = Channel.fromPath( params.aligned_bam )
    aligned_reads_idx = Channel.fromPath( params.aligned_bam + ".bai" )

    medaka_snv_calling( aligned_reads, aligned_reads_idx, genomeref, genomeindex )

}


workflow {

    genomeref = Channel.fromPath( params.genomeref + "/fasta/genome.fa", checkIfExists: true  )
    genomeindex = Channel.fromPath( params.genomeref + "/indexes/minimap2-ont/genome.mmi" )
    // ont_base = Channel.fromPath( params.ont_base_dir, checkIfExists: true )
    reads = Channel.fromPath( params.processed_reads, checkIfExists: true )

    // alignment and covnersion into indexed sorted bam
    minimap_alignment( genomeref, genomeindex, reads )
    sam_to_sorted_bam( minimap_alignment.out.mapped_sam, genomeref, minimap_alignment.out.aligner )
    sv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

    // process_reads( genomeref, ont_base )
    // lra_alignment_sv_calling( process_reads.out.fastq_trimmed, genomeref )
    // minimap_alignment_snv_calling( process_reads.out.fastq_trimmed, genomeref )

}

// -m ${Math.floor( $task.memory / $task.cpus ) }
// ch_reference_fasta.view()
// ch_reference_index.view()