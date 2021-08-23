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

/* 
 * pipeline input parameters 
 */
// params {
params.genomeref           = "/staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chr20/"
// params.genomerefindex      = "${file(params.genomeref).getParent()}/${file(params.genomeref).getBaseName()}_map-ont.mmi"
// params.fastqs              = "/staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/data/20210503_S2_2folddilser_100kto3k_nano-gTag/20210503_1530_MN34250_AGI654_ad1ed051/fastq_pass/barcode06/"
params.ont_base_dir        = "/staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/data/raw/testsample/20210719_ASA_Edin_BA24_14_18/"
params.guppy_gpu           = true
params.min_read_qscore     = 10
params.sampleid            = "ASA_Edin_BA24_14_18_chr20"
params.outdir              = "/staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/data/mapped/${params.sampleid}"
params.tracedir            = "${params.outdir}/pipeline_info"
params.sv_min_support      = 8
params.sv_min_mapq         = 10
params.sv_min_size         = 30
params.megalodon_model     = "res_dna_r941_min_modbases-all-context_v001.cfg"
params.medaka_snp_model    = "r941_prom_sup_snp_g507"
params.guppy_config        = "dna_r9.4.1_450bps_sup_prom.cfg"
// params.rerio_base          = "/staging/leuven/stg_00002/lcb/jdemeul/software/rerio/basecall_models/"
// params.guppy_barcode_kit   = "/staging/leuven/stg_00002/lcb/jdemeul/software/rerio/basecall_models/"
params.with_gpu            = true
params.gpu_devices         = "3"
params.rebasecall          = false
params.processed_reads     = ""
params.aligned_bam         = ""
// }


/* 
* Basecall reads
*/
process basecall_reads {
    label 'process_middle'
    label 'guppy'
    label ( params.with_gpu ? 'with_gpus': null )

    publishDir path: "./results/basecalls/", mode: 'copy'

    input:
    path ont_base
    path genomeref
    path genomrefidx

    output:
    path "basecalls/*fastq.gz", emit: fastqs
    path "basecalls/sequencing_summary.txt"

    script:
    """
    guppy_basecaller \
        -i $ont_base \
        -s basecalls \
        -c ${params.guppy_config} \
        --recursive \
        --device "cuda:${params.gpu_devices}" \
        --align_ref $genomeref \
        --compress_fastq \
        --disable_qscore_filtering
    """
}



/* 
* Process fastq files
*/
process filter_reads {
    label 'process_medium'
    label 'fastp'

    publishDir path: "./results/fastq/", mode: 'copy'

    input:
    path fastqs

    output:
    path "*_trimmed.fastq.gz", emit: fastq_trimmed
    path "fastp*"

    script:
    """
    zcat $fastqs | fastp --stdin \
        --disable_adapter_trimming \
        --average_qual ${params.min_read_qscore} --qualified_quality_phred 0 --unqualified_percent_limit 100 --n_base_limit 50 \
        --length_required 100 \
        --trim_front1 30 --trim_tail1 15 \
        -o ${params.sampleid}_trimmed.fastq.gz \
        --thread $task.cpus 
    """
}



/* 
* index a reference genome with minimap2
*/
process create_minimap_index {
    // tag "$genomeref"
    label 'process_medium'
    label 'minimap'

    // publishDir path: '.', mode: 'copy',
    //     saveAs: { "${file(params.genomerefindex)}" }
    publishDir path: "${file(params.genomeref).getParent() + '/indexes/minimap2-ont/'}", mode: 'copy'
        // saveAs: { "${file(params.genomerefindex)}" }

    // when:
    // !file(params.genomerefindex).exists()
    // genomeref.getExtension() != "mmi"

    input:
    path genomeref // from ch_reference_fasta

    output:
    path "genome.mmi", emit: mmi // into ch_new_reference_index

    script:
    // if ( !file(params.genomerefindex).exists() )
    """
    minimap2 -x map-ont -k 17 -t $task.cpus -d genome.mmi $genomeref
    """
    // else 
    //     """
    //     cp ${params.genomerefindex} index
    //     """

}


/* 
* align reads to the reference with minimap2
*/
process minimap_alignment {
    label 'process_high'
    label 'minimap'

    input:
    path genomeref
    path index // from ch_old_reference_index.mix(ch_new_reference_index)
    path reads // from ch_fastqs

    output:
    path "mapped.sam", emit: mapped_sam
    val "minimap2", emit: aligner

    script:
    def samtools_mem = Math.floor(task.memory.getMega() / task.cpus ) as int
    """
    minimap2 -ax map-ont -t $task.cpus -L --secondary=no --MD --cap-kalloc=500m -K 5g $genomeref ${reads} > mapped.sam
    """
}


/* 
* index a reference genome with LRA
*/
process create_lra_index {
    // tag "$genomeref"
    label 'process_medium'
    label 'lra'

    publishDir path: "${file(params.genomeref).getParent() + '/indexes/lra-ont/'}", mode: 'copy'

    input:
    path genomeref // from ch_reference_fasta

    output:
    path "*.gli", emit: gli // into ch_new_reference_index
    path "*.mmi", emit: mmi // into ch_new_reference_index

    script:
    """
    lra index -ONT $genomeref
    """
}

/* 
* align reads to the reference with LRA
*/
process lra_alignment {
    label 'process_high'
    label 'lra'

    input:
    path genomeref
    // path genomeindex
    path reads // from ch_fastqs
    path gliidx
    path mmiidx

    output:
    path "mapped.sam", emit: mapped_sam
    val "lra", emit: aligner

    script:
    """
    zcat ${reads} | lra align -ONT -p s -t $task.cpus $genomeref /dev/stdin > mapped.sam
    """
}


/* 
* sam to bam conversion using samtools
*/
process sam_to_sorted_bam {
    label 'process_medium'
    label 'samtools'

    publishDir path: "./results/bam/", mode: 'copy',
               saveAs: { item -> item.matches("(.*)minimap2(.*)") ? item : null }

    input:
    path mapped_sam
    path genomeref
    val aligner

    output:
    path "*.bam", emit: sorted_bam
    path "*.bai", emit: bam_index
    path "*stats"

    script:
    def samtools_mem = Math.floor(task.memory.getMega() / task.cpus ) as int
    """
    samtools sort -@ $task.cpus \
        --write-index \
        -o ${params.sampleid}_${aligner}.bam##idx##${params.sampleid}_${aligner}.bam.bai \
        -m ${samtools_mem}M \
        --reference $genomeref \
        -T sorttmp_${params.sampleid}_${aligner} \
        $mapped_sam
    samtools flagstat ${params.sampleid}_${aligner}.bam > ${params.sampleid}_${aligner}.bam.flagstats
    samtools idxstats ${params.sampleid}_${aligner}.bam > ${params.sampleid}_${aligner}.bam.idxstats
    samtools stats ${params.sampleid}_${aligner}.bam > ${params.sampleid}_${aligner}.bam.stats
    """

}


/* 
* SV calling on a bam file using cuteSV
*/
process cutesv_sv_calling {
    label 'process_high'
    label 'cutesv'

    publishDir path: "./results/svs_cutesv/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "*_cuteSV_svs.vcf", emit: sv_calls

    script:
    """
    cuteSV --min_support ${params.sv_min_support} \
        --min_mapq ${params.sv_min_mapq} \
        --report_readid \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --max_cluster_bias_INV 500 \
        --max_cluster_bias_DUP 500 \
        --max_cluster_bias_TRA 50 \
        --diff_ratio_filtering_TRA 0.6 \
        --genotype \
        --min_size ${params.sv_min_size} \
        --sample ${params.sampleid} \
        --threads $task.cpus \
        $sorted_bam \
        $genomeref \
        ${params.sampleid}_cuteSV_svs.vcf \
        `pwd`
    """

}


/* 
* SV calling on a bam file using Sniffles
*/
process sniffles_sv_calling {
    label 'process_high'
    label 'sniffles'

    publishDir path: "./results/svs_sniffles/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index

    output:
    path "*_sniffles_svs.vcf", emit: sv_calls

    script:
    """
    sniffles --minmapping_qual ${params.sv_min_mapq} \
        --min_support ${params.sv_min_support} \
        --cluster \
        --min_length ${params.sv_min_size} \
        --num_reads_report -1 \
        --threads $task.cpus \
        --tmp_file working_dir \
        -m $sorted_bam \
        -v ${params.sampleid}_sniffles_svs.vcf
    """

}


/* 
* SV calling on a bam file using SVIM
*/
process svim_sv_calling {
    label 'process_low'
    label 'svim'

    publishDir path: "./results/svs_svim/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "working_dir/*.vcf", emit: sv_calls
    path "working_dir/*.png"

    script:
    """
    svim alignment \
        --min_mapq ${params.sv_min_mapq} \
        --min_sv_size ${params.sv_min_size} \
        --sample ${params.sampleid} \
        --read_names \
        working_dir \
        $sorted_bam \
        $genomeref
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process svim_sv_filtering {
    label 'process_low'
    label 'bcftools'

    publishDir path: "./results/svs_svim/", mode: 'copy'

    input:
    path svs

    output:
    path "*_Q10.vcf", emit: sv_calls_q10

    script:
    """
    bcftools view -i 'QUAL >= 10' -o ${svs.getSimpleName()}_Q10.vcf $svs
    """
}



/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process deepvariant_snv_calling {
    label 'process_high'
    label 'deepvariant'

    publishDir path: "./results/snv_indel_deepvariant/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "*.vcf.gz", emit: indel_snv_vcf
    path "*.vcf.gz.tbi", emit: indel_snv_vcf_index
    path "logs"
    path "intermediate_files"
    path "*.visual_report.html"

    script:
    """
    run_pepper_margin_deepvariant call_variant \
        -b $sorted_bam \
        -f $genomeref \
        -o . \
        -p ${params.sampleid} \
        -t $task.cpus \
        --ont \
        --phased_output
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process medaka_snv_calling {
    label 'process_high'
    label 'medaka'
    label (${params.with_gpu} ? 'with_gpu': null)

    publishDir path: "./results/snv_indel_medaka/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "medaka_variant", emit: medaka_variant

    script:
    """
    export SINGULARITYENV_CUDA_VISIBLE_DEVICES=${params.gpu_devices}
    medaka_variant \
        -i $sorted_bam \
        -f $genomeref \
        -t $task.cpus \
        -b 150 \
        -n ${params.sampleid} \
        -s ${params.medaka_snp_model} \
        -m ${params.medaka_snp_model.replace("snp", "variant")} \
        -o medaka_variant \
        -p 
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process variant_filtering {
    label 'process_high'
    label 'bcftools'

    publishDir path: "./results/", mode: 'copy'
    // ,
        // saveAs: { item -> 
        //                 if ( item.matches("(.*)lra(.*)" ) {
        //                     "/svs/" + item
        //                 } else if item.matches("(.*)medaka(.*)" ) {
        //                     "/snvs_indel_medaka/" + item
        //                 } else {
        //                     "/snvs_indel/" + item
        //                 }}

    input:
    path variants

    output:
    path "*_PASS.vcf", emit: variants_pass

    script:
    """
    bcftools view -f "PASS" -o ${variants.getSimpleName()}_PASS.vcf $variants
    """
}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process create_personal_genome {
    label 'process_low'
    label 'deepvariant'

    publishDir path: "./results/crossstitch", mode: 'copy'

    input:
    path phased_snps
    path unphased_svs
    path sorted_bam
    path bam_index
    path genomeref
    val karyotype
    val refine

    output:
    path "*.vcf.gz", emit: indel_snv_vcf
    path "*.vcf.gz.tbi", emit: indel_snv_vcf_index
    path "logs"
    path "intermediate_files"
    path "*.visual_report.html"

    script:
    """
    export XDG_CONFIG_HOME="/staging/leuven/stg_00002/lcb/jdemeul/software/"
    module use /staging/leuven/stg_00002/lcb/jdemeul/software/easybuild/modules/all/
    ml --ignore-cache Java/13.0.2 SAMtools/1.9-GCC-6.4.0-2.28 HTSlib/1.9-GCC-6.4.0-2.28 pigz/2.6-GCCcore-6.4.0 Racon/1.4.13-GCCcore-9.3.0
    /staging/leuven/stg_00002/lcb/jdemeul/software/crossstitch/src/crossstitch.sh  \
        /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20.phased_PASS.vcf /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20_LRA_cuteSV_svs_PASS.vcf /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/code/nanowgs/results/ASA_Edin_BA24_14_18_chr20_minimap2.bam /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chr20.fasta ASA_Edin_BA24_14_18_chr20 xy 1
    run_pepper_margin_deepvariant call_variant \
        -b $sorted_bam \
        -f $genomeref \
        -o . \
        -p ${params.sampleid} \
        -t $task.cpus \
        --ont \
        --phased_output
    """
}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process megalodon_modifications {
    label 'process_high'
    label 'megalodon'
    label (${params.with_gpu} ? 'with_gpu': null)

    publishDir path: "./results/", mode: 'copy'

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


/* 
* Run de novo genome assembly using Shasta
*/
process run_shasta_assembly {
    label 'bigmen'
    label 'shasta'

    input:
    path phased_snps

    output:
    path "*.vcf.gz", emit: indel_snv_vcf

    script:
    """
    singularity run \
    -B /staging/leuven/stg_00002/lcb/jdemeul/ \
    -B /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/shasta_assembly/:/output \
    /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-shasta-docker-latest.img 0.7.0 \
    --input /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ASA_Edin_BA24_38_17_trimmed.fastq \
    --conf /staging/leuven/stg_00002/lcb/jdemeul/software/shasta/conf/Nanopore-Sep2020.conf
    """
}


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


workflow lra_alignment_sv_calling {
    take: 
        fastqs
        genomeref
    main:
        // genomeref = Channel.fromPath( params.genomeref + "/fasta/genome.fa", checkIfExists: true  )
        // fastqs = Channel.fromPath( params.fastqs )

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


workflow process_reads {
    take:
        genomeref
        ont_base

    main:
    
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