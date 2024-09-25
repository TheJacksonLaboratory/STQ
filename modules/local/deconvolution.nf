
process XENOME_GENERATE_INDEX {

    tag "${params.deconvolution_indices_name}"
    publishDir "${params.deconvolution_indices_path}", pattern: "${params.deconvolution_indices_name}-*", mode: 'copy', overwrite: false

    input:
    path host_fasta
    path graft_fasta
    val kmer_size
    
    output:
    path("${params.deconvolution_indices_name}-*"), emit: indices_path
    
    script:    
    """
    mkdir tempw
    
    /xenome-1.0.1-r/xenome index \
    --kmer-size ${kmer_size} \
    --prefix ${params.deconvolution_indices_name} \
    --tmp-dir tempw \
    --num-threads ${task.cpus} \
    --host "${host_fasta}" \
    --graft "${graft_fasta}" \
    --verbose \
    --max-memory 20
    """
}

process XENGSORT_GENERATE_INDEX {

    tag "${params.deconvolution_indices_name}"
    publishDir "${params.deconvolution_indices_path}", pattern: "${params.deconvolution_indices_name}-xind*", mode: 'copy', overwrite: false

    input:
    path host_fasta
    path graft_fasta
    val kmer_size
    
    output:
    path("${params.deconvolution_indices_name}-xind*"), emit: indices_path
    
    script:    
    """    
    xengsort index \
    -H "${host_fasta}" \
    -G "${graft_fasta}" \
    -n ${params.xengsort_n} \
    -k ${kmer_size} \
    -W ${task.cpus} \
    --index ${params.deconvolution_indices_name}-xind
    """
}

process DECONVOLUTION_XENOME {

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

process DECONVOLUTION_XENGSORT {

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", pattern: '.command.out', saveAs: { filename -> "xengsort.summary.txt" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastq)
    path(indices_path)
    val(indices_name)
    
    output:
    tuple val(sample_id), file("unsorted_human_{1,2}.fastq"), emit: human
    tuple val(sample_id), file("unsorted_mouse_{1,2}.fastq"), emit: mouse
    path(".command.out"), emit: summary
    
    script:    
    """
    xengsort classify \
    --index "${indices_path}/${params.deconvolution_indices_name}-xind" \
    --fastq ${fastq[1]} \
    --out fastq \
    --threads ${task.cpus} \
    --chunksize 32.0 \
    --compression none

    mv fastq-host.fq unsorted_mouse_2.fastq
    mv fastq-graft.fq unsorted_human_2.fastq
    
    seqtk subseq ${fastq[0]} <(fgrep "@" unsorted_mouse_2.fastq | cut -d ' ' -f1 | cut -d '@' -f2) >> unsorted_mouse_1.fastq
    seqtk subseq ${fastq[0]} <(fgrep "@" unsorted_human_2.fastq | cut -d ' ' -f1 | cut -d '@' -f2) >> unsorted_human_1.fastq
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
