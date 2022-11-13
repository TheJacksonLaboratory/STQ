
process TEST_VISIBLE {

    label 'process_hovernet'

    script:
    """
    echo Allocated device: $CUDA_VISIBLE_DEVICES
    """ 
}