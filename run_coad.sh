samplesheet="./assets/samplesheet_coad_part3.csv"
workdir="/fastscratch/domans/work_coad_part3"
outdir="./results_coad_part3"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
