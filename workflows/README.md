
## Routes of analysis

<p>
    <img src="../docs/route-map.png" width="800"/>
</p>

All three routes of analysis are implemented as Nextflow DSL2 workflows and use the same style [samplesheet](../README.md#samplesheet). The name of the workflow (one of "two_references", "one_reference", and "arbitrary_grid") is specified when the pipeline is invoked:

      nextflow run main.nf [...] --workflow="two_references"

JAX users can edit the file [submit.sb](../submit.sb) to specify the workflow.
