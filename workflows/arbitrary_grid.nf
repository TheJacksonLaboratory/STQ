
include { IMG } from '../subworkflows/imaging'

workflow ARB {

    take:
        samples

    main:
        IMG ( samples
              .join(samples.map{[it[0], []]}))

}
