samplesheet="./assets/samplesheet_3_ADD.csv"
workdir="/fastscratch/domans/work_ADD"
outdir="./results_3_AD"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
