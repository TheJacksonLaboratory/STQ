#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TEST_VISIBLE {

    label 'process_hovernet'

    script:
    """
    echo Allocated device: $CUDA_VISIBLE_DEVICES
    """ 
}


process TEST_RETRIES {
        
    maxRetries params.maxRetries
    errorStrategy { task.attempt > params.maxRetries ? 'ignore' : 'retry' }
    
    debug true
    
    script: 
    println task.attempt + ", " +  task.maxRetries + ", " +  task.errorStrategy
       
    """    
    echo ${task.attempt} of ${task.maxRetries}, errorStrategy: ${task.errorStrategy}
    
    I-am-groot
    """ 
}


workflow {

    TEST_RETRIES()

}