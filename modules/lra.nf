
/* 
* Index a reference genome with LRA
*/
process create_lra_index {
    // tag "$genomeref"
    label 'process_medium'
    label 'lra'

    publishDir path: "${file(params.genomeref).getParent() + '/indexes/lra-ont/'}", mode: 'copy'

    input:
    path genomeref

    output:
    path "*.gli", emit: gli
    path "*.mmi", emit: mmi

    script:
    """
    lra index -ONT $genomeref
    """
}

/* 
* Align reads to a reference genome with LRA
*/
process lra_alignment {
    label 'process_high'
    label 'lra'

    input:
    path genomeref
    path reads
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
