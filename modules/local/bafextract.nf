
process GET_REFERENCE_PILEUP {

    tag "$reference"
    label "bafextract"

    input:
    path(reference)
    
    output:
    tuple file("sizes.list"), file("pileup/*")
    
    script:    
    """
    mkdir pileup

    /BAFExtract/bin/BAFExtract -preprocess_FASTA ${reference}/fasta/genome.fa pileup
    
    cut -f1,2 ${reference}/fasta/genome.fa.fai > sizes.list 
    """
}


process GET_PILEUP_OF_BAM {

    tag "$sample_id"
    label "bafextract"

    input:
    tuple val(sample_id), path(bam)
    tuple path(sizes), path(pileup)
    
    output:
    tuple val(sample_id), file("bam_pileup/*")
    
    script:    
    """
    mkdir bam_pileup

    if [[ ${bam[0]} == *".bam.bai"* ]]; then
        bbam=${bam[1]}
    else
        bbam=${bam[0]}
    fi

    samtools view \${bbam} | /BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin ${sizes} bam_pileup ${params.bafextract_minimum_mapping_quality} ${params.bafextract_minimum_base_quality}
    """
}


process GET_SNV_FROM_PILEUP {

    tag "$sample_id"
    label "bafextract"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(bam_pileup)
    tuple path(sizes), path(pileup)
    val(species)
    
    output:
    tuple val(sample_id), file("${species}/extracted.baf")
    
    script:    
    """
    mkdir ref_pileup
    
    for f in ${pileup}
    do
        if [[ \${f} == *" "* ]]
        then
            newf="`cut -d ' ' -f 1 <<< "\$f"`.bin"
            cp "\${f}" "ref_pileup/\${newf}"
        else
            cp "\${f}" "ref_pileup/\${f}"
        fi
    done
      
    mkdir bam_pileup
    cp ${bam_pileup} bam_pileup
    
    mkdir "${species}"
    
    /BAFExtract/bin/BAFExtract -get_SNVs_per_pileup ${sizes} bam_pileup ref_pileup ${params.bafextract_min_coverage_per_SNV} ${params.bafextract_min_MAF_covg_per_SNV} ${params.bafextract_min_MAF} "${species}/extracted.baf"
    
    rm -R ref_pileup
    rm -R bam_pileup
    """
}
