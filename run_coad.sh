samplesheet="./assets/samplesheet_tcga_coad_10000_p3.csv"
workdir="/fastscratch/domans/work_coad_10000"
outdir="./results_coad_10000"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
