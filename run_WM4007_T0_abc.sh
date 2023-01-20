samplesheet="./assets/samplesheet_WM4007_3_T0_abc.csv"
workdir="/fastscratch/domans/work_ADD"
outdir="./results_3_AD"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
