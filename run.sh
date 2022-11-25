samplesheet="./assets/samplesheet_WM4237_ADD.csv"
workdir="/fastscratch/domans/work_WM4237_ADD"
outdir="./result"

sbatch job_sumner.sb ${samplesheet} ${workdir} ${outdir}
