
workflow="arbitrary_grid" ### "two_references" "one_reference" "arbitrary_grid" "deconvolution_indices"
samplesheet="./assets/samplesheet-kidney-wedge-test.csv"
workdir="/flashscratch/domans/work"
outdir="../results-test"
binddir="/projects/"

#----------------------------------------------------------------------------------------------------

./check.sh

SLURM_SUBMIT_DIR=`pwd`

sbatch \
submit.sb $workflow "$samplesheet" "$workdir" "$outdir" "$binddir"
