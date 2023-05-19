samplesheet="./assets/samplesheet_test.csv"
workdir="/fastscratch/domans/work_test"
outdir="./results_test"

job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
#sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
