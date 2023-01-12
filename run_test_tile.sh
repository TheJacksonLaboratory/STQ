samplesheet="./assets/samplesheet_test_tile.csv"
workdir="/fastscratch/domans/work_test_tile"
outdir="./results_test_tile"

#sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
