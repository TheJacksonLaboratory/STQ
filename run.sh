samplesheet="./assets/samplesheet.csv"
workdir="/fastscratch/domans/work_PDXq_test"
outdir="./results_test"

sbatch job_sumner.sb ${samplesheet} ${workdir} ${outdir}
