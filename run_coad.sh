samplesheet="./assets/samplesheet_tcga_coad_test.csv"
workdir="/fastscratch/domans/work_coad_test"
outdir="./results_coad_test"

sbatch job_sumner.sb "${samplesheet}" "${workdir}" "${outdir}"
