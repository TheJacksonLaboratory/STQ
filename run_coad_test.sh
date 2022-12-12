samplesheet="./assets/samplesheet_test_coad.csv"
workdir="/fastscratch/domans/work_test_coad"
outdir="./results_test_coad"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
