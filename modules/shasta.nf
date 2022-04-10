
/* 
* De novo genome assembly using Shasta
*/
process run_shasta_assembly {
    label 'shasta'
    label ( workflow.profile.contains('qsub') ? 'bigmem': 'cpu_high' )
    label ( workflow.profile.contains('qsub') ? null: 'mem_high' )
    label ( workflow.profile.contains('qsub') ? null: 'time_mid' )

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path fastq
    // path config

    output:
    path "shasta_assembly/Assembly.fasta", emit: assembly
    path "shasta_assembly"

    script:
    def localproc = ( workflow.profile.contains('qsub') ? 0: task.cpus )
    if( params.shasta_minreadlength )
        """
        if [[ $fastq == *.gz ]]; then 
            gunzip -c $fastq > uncompressed_reads.fq
        else 
            mv $fastq uncompressed_reads.fq
        fi

        shasta \
            --config /shastaconf/conf/${params.shasta_config} \
            --input uncompressed_reads.fq \
            --assemblyDirectory ./shasta_assembly \
            --threads ${localproc} \
            --Reads.minReadLength ${params.shasta_minreadlength}
        
        if [ -f ./shasta_assembly/Assembly-Haploid.fasta ]
        then mv ./shasta_assembly/Assembly-Haploid.fasta ./shasta_assembly/Assembly.fasta
        fi 
        """
    else
        """
        if [[ $fastq == *.gz ]]; then 
            gunzip -c $fastq > uncompressed_reads.fq
        else 
            mv $fastq uncompressed_reads.fq
        fi

        shasta \
            --config /shastaconf/conf/${params.shasta_config} \
            --input uncompressed_reads.fq \
            --assemblyDirectory ./shasta_assembly \
            --threads ${localproc}
        
        if [ -f ./shasta_assembly/Assembly-Haploid.fasta ]
        then mv ./shasta_assembly/Assembly-Haploid.fasta ./shasta_assembly/Assembly.fasta
        fi 
        """
}


/* 
* De novo haploid genome assembly using Shasta
* DEPRECATED
*/
process run_shasta_assembly_haploid {
    label 'shasta'
    label ( workflow.profile.contains('qsub') ? 'bigmem': 'cpu_high' )
    label ( workflow.profile.contains('qsub') ? null: 'mem_high' )
    label ( workflow.profile.contains('qsub') ? null: 'time_mid' )

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process}/", mode: 'copy'

    input:
    path hapreads
    val haplotype

    output:
    path "${haplotype}/Assembly.fasta", emit: assembly
    path "${haplotype}"

    script:
    def localproc = ( workflow.profile.contains('qsub') ? 0: task.cpus )
    if( params.shasta_minreadlength )
        """
        shasta \
            --config /shastaconf/conf/${params.shasta_config_haploid} \
            --input $hapreads \
            --assemblyDirectory $haplotype \
            --threads ${localproc} \
            --Reads.minReadLength ${params.shasta_minreadlength}
        """
    else
        """
        shasta \
            --config /shastaconf/conf/${params.shasta_config_haploid} \
            --input $hapreads \
            --assemblyDirectory $haplotype \
            --threads ${localproc}
        """
}



    // singularity run -B /staging/leuven/stg_00002/lcb/jdemeul/ /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-mummer-4.0.0rc1.img nucmer -t 16 --maxmatch -l 100 -c 1000 -p /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/fasta/chm13_v1.1_chrY_KI270740_EBV.fasta /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/Assembly.fasta
    // python DotPrep.py --delta /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ASA_Edin_BA24_38_17_denovo_vs_CHM13v1.1.delta --out /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ASA_Edin_BA24_38_17_denovo_vs_CHM13v1.1.DotPrep.out
 
