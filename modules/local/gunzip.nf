
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


process GUNZIP_SEPARATE {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), file("*_R1_*.fastq"), emit: R1
    tuple val(sample_id), file("*_R2_*.fastq"), emit: R2
    
    script:
    """   
    gzip -d -k ${fastq}
    """
}


process GUNZIP_FASTA {

    input:
    path(fasta)
    
    output:
    file("*{.fa,.fna}")
    
    script:
    """   
    gzip -d -k ${fasta}
    """
}
