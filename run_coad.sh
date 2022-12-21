samplesheet="./assets/samplesheet_test_coad.csv"
workdir="/fastscratch/domans/work_coad_test"
outdir="./results_coad_test"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
