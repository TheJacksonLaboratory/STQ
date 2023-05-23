
process LOAD_SAMPLE_INFO {

    tag "$sample_id"

    input:
    tuple val(sample_id), val(meta), path(fastq), path(image)
    
    output:
    tuple val(sample_id), file(image), emit: image
    tuple val(sample_id), file("${fastq}/*R{1,2}*.fastq*"), emit: fastq

    script:    
    """
    """
}
