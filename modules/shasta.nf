
/* 
* De novo genome assembly using Shasta
*/
process run_shasta_assembly {
    label 'shasta'
    label ( workflow.profile.contains('qsub') ? 'bigmem': 'cpu_high' )
    label ( workflow.profile.contains('qsub') ? null: 'mem_high' )
    label ( workflow.profile.contains('qsub') ? null: 'time_mid' )

    publishDir path: "${params.outdir}/results/", mode: 'copy'

    input:
    path fastq
    path config

    output:
    path "shasta_assembly/Assembly.fasta", emit: assembly
    path "shasta_assembly/*"

    script:
    """
    if [[ $fastq == *.gz ]]; then 
        gunzip -c $fastq > uncompressed_reads.fq
    else 
        mv $fastq uncompressed_reads.fq
    fi

    /opt/shasta-Linux-0.7.0 \
        --input uncompressed_reads.fq \
        --conf $config \
        --assemblyDirectory ./shasta_assembly
    """
}



    // singularity run -B /staging/leuven/stg_00002/lcb/jdemeul/ /staging/leuven/stg_00002/lcb/jdemeul/software/singularity_images/zeunas-mummer-4.0.0rc1.img nucmer -t 16 --maxmatch -l 100 -c 1000 -p /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ /staging/leuven/stg_00002/lcb/jdemeul/reference/chm13_v1.1_chrY_KI270740_EBV/fasta/chm13_v1.1_chrY_KI270740_EBV.fasta /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/Assembly.fasta
    // python DotPrep.py --delta /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ASA_Edin_BA24_38_17_denovo_vs_CHM13v1.1.delta --out /staging/leuven/stg_00002/lcb/jdemeul/projects/2021_ASAP/results/ASA_Edin_BA24_38_17/results/fastq/ShastaRun/mummer_ASA_Edin_BA24_38_17_denovo_vs_CHM13v1.1.DotPrep.out
 
