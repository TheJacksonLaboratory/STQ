
process CONVERT_TO_PYRAMIDAL_OME {

    tag "$sample_id"
    label 'process_ome'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'image.ome.tiff', mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("image.ome.tiff")

    script:    
    """
    export BF_MAX_MEM=24g
    export _JAVA_OPTIONS="-Xmx24g"
    bfconvert -version

    bfconvert -noflat -bigtiff -overwrite \
    -pyramid-resolutions 3 -pyramid-scale 4 -tilex ${params.tiled_tiff_tile_size} -tiley ${params.tiled_tiff_tile_size} \
    -compression ${params.compression} "${image}" image.ome.tiff || \
    bfconvert -noflat -bigtiff -overwrite \
    -pyramid-resolutions 2 -pyramid-scale 4 -tilex ${params.tiled_tiff_tile_size} -tiley ${params.tiled_tiff_tile_size} \
    -compression ${params.compression} "${image}" image.ome.tiff || \
    bfconvert -noflat -bigtiff -overwrite \
    -pyramid-resolutions 1 -pyramid-scale 4 -tilex ${params.tiled_tiff_tile_size} -tiley ${params.tiled_tiff_tile_size} \
    -compression ${params.compression} "${image}" image.ome.tiff
    """
}


process EXTRACT_IMAGE_METADATA {

    tag "$sample_id"
    label 'process_ome'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'metadata.ome.xml', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 4.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 0.GB }

    input:
    tuple val(sample_id), path(fileslide), path(roifile), val(mpp), val(size)
    
    output:
    tuple val(sample_id), file("metadata.ome.xml")

    script:    
    """
    showinf -omexml-only -nopix ${fileslide} >> "metadata.ome.xml"
    """
}