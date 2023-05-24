
process XENOME_GENERATE_INDEX {

    tag "${params.xenome_indices_name}"
    publishDir "${params.xenome_indices_path}", pattern: "${params.xenome_indices_name}-*", mode: 'copy', overwrite: false

    input:
    path host_fasta
    path graft_fasta
    val kmer_size
    
    output:
    path("${params.xenome_indices_name}-*"), emit: indices_path
    
    script:    
    """
    /xenome-1.0.1-r/xenome index \
    --kmer-size ${kmer_size} \
    --prefix ${params.xenome_indices_name} \
    --tmp-dir tempinwork
    --num-threads ${task.cpus} \
    --host "${host_fasta}" \
    --graft "${graft_fasta}" \
    --verbose
    """
}

process DECONVOLUTION {

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", pattern: '.command.out', saveAs: { filename -> "xenome.summary.txt" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq)
    path(indices_path)
    val(indices_name)
    
    output:
    tuple val(sample_id), file("categorized/unsorted_human_{1,2}.fastq"), emit: human
    tuple val(sample_id), file("categorized/unsorted_mouse_{1,2}.fastq"), emit: mouse
    path(".command.out"), emit: summary
    
    script:    
    """
    mkdir categorized
    mkdir tmp

    /xenome-1.0.1-r/xenome classify \
    -T ${task.cpus} \
    -i ${fastq[0]} \
    -i ${fastq[1]} \
    --pairs \
    -P ${indices_path}/${indices_name} \
    --graft-name human \
    --host-name mouse \
    --output-filename-prefix categorized/unsorted \
    --tmp-dir tmp \
    --verbose
    """
}

process SORT_FASTQ {

    tag "$sample_id"
    label "low_process"

    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), file("sorted_{1,2}.fastq")
    
    script:    
    """
    fastq-sort --id ${fastq[0]} > "sorted_1.fastq"
    fastq-sort --id ${fastq[1]} > "sorted_2.fastq"
    """
}
