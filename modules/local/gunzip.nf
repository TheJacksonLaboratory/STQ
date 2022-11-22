
process GUNZIP {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), file("*R{1,2}*.fastq")
    
    script:
    """   
    gzip -d -k ${fastq}
    """
}
