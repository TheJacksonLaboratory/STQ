samplesheet="./assets/samplesheet_WM4007_WM4237_AD_v3.csv"
workdir="/fastscratch/domans/work_WM4007_WM4237_AD"
outdir="./results_WM4007_WM4237_AD"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
