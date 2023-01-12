samplesheet="./assets/samplesheet_demo.csv"
workdir="/fastscratch/domans/work_demo"
outdir="./results_demo"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
