samplesheet="./assets/samplesheet_test_aperio.csv"
workdir="/fastscratch/domans/work_test_aperio"
outdir="./results_test_aperio"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
