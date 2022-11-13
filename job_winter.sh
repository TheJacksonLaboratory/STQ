stage=$1
samplesheet=$2
workdir=$3

module use --append /projects/omics_share/meta/modules
module load nextflow

nextflow run main.nf \
-w ${workdir} \
-profile slurm,singularity \
-resume \
--input=${samplesheet} \
--stage=${stage}
