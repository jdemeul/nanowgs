
/* 
* SNV and indel calling on aligned reads using PEPPER-Margin-DeepVariant
*/
process deepvariant_snv_calling {
    // label 'process_high'
    label 'deepvariant'
    label ( params.with_gpu ? 'with_gpus': null )

    publishDir path: "${params.outdir}/results/snv_indel_deepvariant/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path genomeref
    // path genomerefidx

    output:
    path "*.vcf.gz", emit: indel_snv_vcf
    path "*.vcf.gz.tbi", emit: indel_snv_vcf_index
    path "logs"
    path "intermediate_files"
    path "*.visual_report.html"

    script:
    if ( params.with_gpu ) 
        """
        export CUDA_VISIBLE_DEVICES=${params.gpu_devices}
        run_pepper_margin_deepvariant call_variant \
            -b $sorted_bam \
            -f $genomeref \
            -o . \
            -p ${params.sampleid} \
            -s ${params.sampleid} \
            -t $task.cpus \
            --ont \
            -g \
            --phased_output
        """
    else 
        """
        run_pepper_margin_deepvariant call_variant \
            -b $sorted_bam \
            -f $genomeref \
            -o . \
            -p ${params.sampleid} \
            -s ${params.sampleid} \
            -t $task.cpus \
            --ont \
            --phased_output
        """
}




/* 
* Assembly polishing using PEPPER
* NOT IMPLEMENTED YET
*/
process pepper_assembly_polishing {
    label 'process_high'
    label 'pepper'
    label ( params.with_gpu ? 'with_gpus': null )

    publishDir path: "${params.outdir}/results/pepper_polished_assembly/", mode: 'copy'

    input:
    path sorted_bam
    path bam_index
    path assemblyref

    output:
    path "*.vcf.gz", emit: indel_snv_vcf

    script:
    """
    run_pepper_margin_deepvariant polish_assembly \
        -b $sorted_bam \
        -f $assemblyref \
        -o . \
        -t $task.cpus \
        -p ${params.sampleid} \
        -g \
        --ont
    """
}


// singularity run --nv -B /staging/leuven/stg_00002/lcb/ \
//     -B /scratch/ \
//     -B /local_scratch/ \
//     /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/kishwars-pepper_deepvariant-r0.4.img \
//     run_pepper_margin_deepvariant polish_assembly \
//     -b /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/PEPPER-Polishing/ASA_Edin_BA24_38_17_TrimmedReads_ShastaAssembly.aln.bam \
//     -f /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/Assembly.fasta \
//     -o /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/PEPPER-Polishing/out \
//     -t 16 \
//     -p ASA_Edin_BA24_38_17 \
//     -g \
//     --ont

//     # this generates 2 VCFs, one per haplotype
//     HAP1_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP1.vcf.gz
//     HAP2_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP2.vcf.gz

//     POLISHED_ASM_HAP1=HG002_Shasta_run1.PMDV.HAP1.fasta
//     POLISHED_ASM_HAP2=HG002_Shasta_run1.PMDV.HAP2.fasta

//     # Apply the VCF to the assembly
//     singularity exec --bind /usr/lib/locale/ \
//     pepper_deepvariant_r0.4.sif \
//     bcftools consensus \
//     -f "${INPUT_DIR}/${ASM}" \
//     -H 2 \
//     -s "${SAMPLE_NAME}" \
//     -o "${OUTPUT_DIR}/${POLISHED_ASM_HAP1}" \
//     "${OUTPUT_DIR}/${HAP1_VCF}"

//     singularity exec --bind /usr/lib/locale/ \
//     pepper_deepvariant_r0.4.sif \
//     bcftools consensus \
//     -f "${INPUT_DIR}/${ASM}" \
//     -H 2 \
//     -s "${SAMPLE_NAME}" \
//     -o "${OUTPUT_DIR}/${POLISHED_ASM_HAP2}" \
//     "${OUTPUT_DIR}/${HAP2_VCF}"
//     """