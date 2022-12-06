samplesheet="./assets/samplesheet_WM4237_ADD.csv"
workdir="/fastscratch/domans/work_WM4237_ADD"
outdir="./results_WM4237_AD"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
