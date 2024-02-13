
workflow="two_references" ### "two_references" "one_reference" "arbitrary_grid" "deconvolution_indices"
samplesheet="./assets/samplesheet_demo.csv"
workdir="./work"
outdir="./results"
binddir="/projects/"

#----------------------------------------------------------------------------------------------------

SLURM_SUBMIT_DIR=`pwd`

sbatch submit.sb $workflow "$samplesheet" "$workdir" "$outdir" "$binddir"
