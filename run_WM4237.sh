samplesheet="./assets/samplesheet_WM4237.csv"
workdir="/fastscratch/domans/work_PDXq_WM4237"
outdir="./results_WM4237"

sbatch job_sumner.sb ${samplesheet} ${workdir} ${outdir}
