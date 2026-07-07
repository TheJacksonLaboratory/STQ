
workflow="cell_mosaic" ### "two_references" "one_reference" "arbitrary_grid" "deconvolution_indices"
samplesheet="./assets/samplesheetTest1.csv"
workdir="./work"
outdir="./results"
binddir="/projects/"

#----------------------------------------------------------------------------------------------------

./check.sh

SLURM_SUBMIT_DIR=`pwd`

# sbatch \
submit.sb $workflow "$samplesheet" "$workdir" "$outdir" "$binddir"
