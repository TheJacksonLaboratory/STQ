
import groovy.json.JsonOutput

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


process EXPORT_PARAMETERS {

    publishDir "${params.tracedir}", pattern: '{*.json}', mode: 'copy', overwrite: params.overwrite_files_on_publish

    output:
    path 'parameters.json'

    script:
    "echo '${JsonOutput.toJson(params)}' > parameters.json"
}


process EXPORT_SAMPLEINFO {

    publishDir "${params.outdir}/${sample_id}", pattern: '{*.json}', mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), val(meta)
    
    output:
    path 'info.json'

    script:
    "echo '${JsonOutput.toJson(meta)}' > info.json"
}
