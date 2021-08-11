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
params.genomeref           = "/staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chr20.fasta"
params.genomerefindex      = "${file(params.genomeref).getParent()}/${file(params.genomeref).getBaseName()}_map-ont.mmi"
// params.fastqs              = "/staging/leuven/stg_00002/lcb/jdemeul/projects/2020_fiberseq/data/20210503_S2_2folddilser_100kto3k_nano-gTag/20210503_1530_MN34250_AGI654_ad1ed051/fastq_pass/barcode06/"
params.fastqs              = "/staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/data/raw/testsample/fastq_pass/"
params.guppy_gpu           = true
params.sampleid            = "ASA_Edin_BA24_14_18_chr20"
params.outdir              = "/staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/data/mapped/${params.sampleid}"
params.tracedir            = "${params.outdir}/pipeline_info"
params.cutesv_min_support  = 8
params.megalodon_model     = "res_dna_r941_min_modbases-all-context_v001.cfg"
params.medaka_snp_model    = "r941_prom_sup_snp_g507"
// }

/* 
* index a reference genome with minimap2
*/
process create_minimap_index {
    // tag "$genomeref"
    label 'process_medium'
    label 'minimap'

    // publishDir path: '.', mode: 'copy',
    //     saveAs: { "${file(params.genomerefindex)}" }
    publishDir path: "${file(params.genomerefindex).getParent()}", mode: 'copy'
        // saveAs: { "${file(params.genomerefindex)}" }

    // when:
    // !file(params.genomerefindex).exists()
    // genomeref.getExtension() != "mmi"

    input:
    path genomeref // from ch_reference_fasta

    output:
    path "*.mmi", emit: refindex // into ch_new_reference_index

    script:
    // if ( !file(params.genomerefindex).exists() )
    """
    minimap2 -ax map-ont -t $task.cpus -d ${genomeref.getBaseName()}_map-ont.mmi $genomeref
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
    path index // from ch_old_reference_index.mix(ch_new_reference_index)
    path reads // from ch_fastqs

    output:
    path "mapped.sam", emit: mapped_sam
    val "minimap2", emit: aligner

    script:
    """
    minimap2 -ax map-ont -t $task.cpus -L --secondary=no $index ${reads}/* > mapped.sam
    """
}


/* 
* index a reference genome with LRA
*/
process create_lra_index {
    // tag "$genomeref"
    label 'process_medium'
    label 'lra'

    // publishDir path: '.', mode: 'copy',
    //     saveAs: { "${file(params.genomerefindex)}" }
    publishDir path: "${file(params.genomerefindex).getParent()}", mode: 'copy'
        // saveAs: { "${file(params.genomerefindex)}" }

    input:
    path genomeref // from ch_reference_fasta

    output:
    path "*.gli" // into ch_new_reference_index
    path "*.mmi" // into ch_new_reference_index

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

    output:
    path "mapped.sam", emit: mapped_sam
    val "lra", emit: aligner

    script:
    """
    zcat ${reads}/*.fastq.gz | lra align -ONT -p s -t $task.cpus $genomeref /dev/stdin > mapped.sam
    """
}


/* 
* sam to bam conversion using samtools
*/
process sam_to_sorted_bam {
    label 'process_medium'
    label 'samtools'

    publishDir path: "./results/", mode: 'copy',
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
        $mapped_sam
    samtools flagstat ${params.sampleid}_${aligner}.bam > ${params.sampleid}_${aligner}.bam.flagstat
    samtools idxstats ${params.sampleid}_${aligner}.bam > ${params.sampleid}_${aligner}.bam.idxstats
    samtools stats ${params.sampleid}_${aligner}.bam > ${params.sampleid}_${aligner}.bam.stats
    """

}


/* 
* SV calling on the LRA bam using cuteSV
*/
process cutesv_sv_calling {
    label 'process_high'
    label 'cutesv'

    publishDir path: "./results/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "*.vcf", emit: sv_calls

    script:
    """
    cuteSV --min_support ${params.cutesv_min_support} \
        --report_readid \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DE 0.3 \
        --genotype \
        --sample ${params.sampleid} \
        --threads $task.cpus \
        $sorted_bam \
        $genomeref \
        ${params.sampleid}_LRA_cuteSV_svs.vcf \
        `pwd`
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process deepvariant_snv_calling {
    label 'process_high'
    label 'deepvariant'

    publishDir path: "./results/", mode: 'copy'

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
        --phased_outputs
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process medaka_snv_calling {
    label 'process_high'
    label 'medaka'

    publishDir path: "./results/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref

    output:
    path "medaka_variant", emit: medaka_variant

    script:
    """
    export CUDA_VISIBLE_DEVICES=2
    medaka_variant \
        -i $sorted_bam \
        -f $genomeref \
        -t $task.cpus \
        -s ${params.medaka_snp_model}
        -m ${params.medaka_snp_model.replace("snp", "variant")}
        -b 150 \
        -o medaka_variant \
        -p 
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process variant_filtering {
    label 'process_high'
    label 'deepvariant'

    publishDir path: "./results/", mode: 'copy'

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
        --phased_outputs
    """

}


/* 
* SNV calling on the minimap2 aligned bam using PEPPER-Margin-DeepVariant
*/
process megalodon_modifications {
    label 'process_high'
    label 'megalodon'

    publishDir path: "./results/", mode: 'copy'

    input:
    path reads
    path variants
    path genomeref

    output:
    path "mod_mappings_sort.bam", emit: mod_basecalls

    script:
    """
    megalodon \
        --devices "cuda:4" \
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



// }

workflow minimap_alignment_snv_calling {
    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    fastqs = Channel.fromPath( params.fastqs )
    
    // genome indexing
    if ( !file(params.genomerefindex).exists() ) {
        create_minimap_index( genomeref )
        genomeindex = create_minimap_index.out.refindex
    } else {
        genomeindex = Channel.fromPath( params.genomerefindex )
    }

    // alignment and covnersion into indexed sorted bam
    minimap_alignment( genomeindex, fastqs )
    sam_to_sorted_bam( minimap_alignment.out.mapped_sam, genomeref, minimap_alignment.out.aligner )

    // SNV calling using PEPPER-margin-DeepVariant
    deepvariant_snv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

}

workflow lra_alignment_sv_calling {
    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true  )
    fastqs = Channel.fromPath( params.fastqs )

    // genome indexing
    if ( !file("${params.genomeref}.gli").exists() ) {
        create_lra_index( genomeref )
    } 

    // alignment and covnersion into indexed sorted bam
    lra_alignment( genomeref, fastqs )
    sam_to_sorted_bam( lra_alignment.out.mapped_sam, genomeref, lra_alignment.out.aligner )

    // SV calling using cuteSV
    cutesv_sv_calling( sam_to_sorted_bam.out.sorted_bam, sam_to_sorted_bam.out.bam_index, genomeref )

}

workflow {
    
    lra_alignment_sv_calling()
    minimap_alignment_snv_calling()

}

// -m ${Math.floor( $task.memory / $task.cpus ) }
// ch_reference_fasta.view()
// ch_reference_index.view()