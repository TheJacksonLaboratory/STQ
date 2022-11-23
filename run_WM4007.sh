samplesheet="./assets/samplesheet_WM4007.csv"
workdir="/fastscratch/domans/work_PDXq_WM4007"
outdir="./results_WM4007"

sbatch job_sumner.sb ${samplesheet} ${workdir} ${outdir}
