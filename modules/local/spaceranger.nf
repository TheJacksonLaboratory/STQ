
process SPACERANGER {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq), path(image)
    path(reference)
    
    output:
    tuple val(sample_id), file("sample/outs/spatial/*"), emit: spatial
    tuple val(sample_id), file("sample/outs/raw_feature_bc_matrix/*"), emit: mtx
    tuple val(sample_id), file("sample/outs/possorted_genome_bam.bam*"), emit: bam
    tuple val(sample_id), file("sample/outs/*summary*"), emit: metrics
    
    script:
    String mem = task.memory
    String memgb = mem.split(" ")[0]    
    """   
    tempfastqdir="temp"
    mkdir \${tempfastqdir}
    cp ${fastq[0]} \${tempfastqdir}/sample_S1_L001_R1_001.fastq
    cp ${fastq[1]} \${tempfastqdir}/sample_S1_L001_R2_001.fastq
    
    spaceranger count \
    --id=sample \
    --sample=sample \
    --fastqs="\${tempfastqdir}" \
    --image=${image} \
    --transcriptome=${reference} \
    --unknown-slide \
    --localcores=${task.cpus} \
    --localmem=${memgb}
    
    rm -R \${tempfastqdir}
    """
}

process RETURN_SPACERANGER_ALIGNMENT {

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), file("mouse/*"), file("human/*"), file("spatial/*")
    
    output:
    tuple val(sample_id), path("spatial/*", includeInputs: true), emit: spatial
    tuple val(sample_id), path("mouse/*", includeInputs: true), path("human/*", includeInputs: true), emit: metrics
    
    script:
    """
    """
}
