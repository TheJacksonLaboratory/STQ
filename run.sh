samplesheet="./assets/samplesheet_WM4237_ADD.csv"
workdir="/fastscratch/domans/work_WM4237_ADD"

sbatch job_sumner.sb ${samplesheet} ${workdir}
