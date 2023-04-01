samplesheet="./assets/samplesheet_tcga_coad_10000_p1.csv"
workdir="/flashscratch/domans/work_coad_10000_p1"
outdir="./results_coad_v3_10000_p1"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
