
process CELLSORT_BAM {

    tag "$sample_id"
    label "samtools"

    input:
    tuple val(sample_id), path("bam/*")
    
    output:
    tuple val(sample_id), file("cellsorted_possorted_genome_bam.bam")
    
    script:    
    """   
    samtools sort bam/possorted_genome_bam.bam -o cellsorted_possorted_genome_bam.bam -t CB -O BAM -@ ${task.cpus}
    """
}


process SPLICING_QUANTIFICATION {

    tag "$sample_id"
    label "splicing_quantification"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path("sample/outs/*"), path("sample/outs/*"), path("sample/outs/filtered_feature_bc_matrix/*")
    path(reference)
    val(species)
    
    output:
    tuple val(sample_id), file("${species}/velocyto.loom")
    
    script:    
    """
    velocyto run10x sample ${reference}/genes/genes.gtf
    
    mkdir ${species}
    cp "sample/velocyto/sample.loom" "${species}/velocyto.loom"
    """
}
