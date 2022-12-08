samplesheet="./assets/samplesheet_WM4007_ADD.csv"
workdir="/fastscratch/domans/work_WM4007_ADD"
outdir="./results_WM4007_AD"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
