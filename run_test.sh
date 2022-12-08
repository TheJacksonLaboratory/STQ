samplesheet="./assets/samplesheet_test.csv"
workdir="/fastscratch/domans/work_test"
outdir="./results_test"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
