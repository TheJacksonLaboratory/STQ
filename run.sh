
workflow="arbitrary_grid" ### "two_references" "one_reference" "arbitrary_grid" "deconvolution_indices"
samplesheet="./assets/samplesheet_test.csv"
workdir="./work"
outdir="./results"
binddir="/projects/"

#----------------------------------------------------------------------------------------------------

SLURM_SUBMIT_DIR=`pwd`

sbatch submit.sb $workflow "$samplesheet" "$workdir" "$outdir" "$binddir"