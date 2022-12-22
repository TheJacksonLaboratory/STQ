samplesheet="./assets/samplesheet_tcga_coad.csv"
workdir="/fastscratch/domans/work_coad"
outdir="./results_coad"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
