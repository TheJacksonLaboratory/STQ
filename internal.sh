#!/usr/bin/bash

stage=$1
samplesheet=$2
workdir=$3

if sacctmgr show clusters | grep -q 'winter';
then
    CLUSTER_NAME="winter"
else
    CLUSTER_NAME="sumner"
fi

if [ $stage == 1 ];
then
    required="sumner"
    if [ $CLUSTER_NAME == $required ];
    then
        echo "Submitting stage ${stage}"
        sbatch submit_${CLUSTER_NAME}.sb ${stage} ${samplesheet} ${workdir}
    else
        echo "ERORR: To run stage ${stage} ssh to $USER@login.${required}.jax.org"
    fi
fi

if [ $stage == 2 ];
then
    required="winter"
    if [ $CLUSTER_NAME == $required ];
    then
        echo "Submitting stage ${stage}"
        sbatch submit_${CLUSTER_NAME}.sb ${stage} ${samplesheet} ${workdir}
    else
        echo "ERORR: To run stage ${stage} ssh to $USER@login.${required}.jax.org"
    fi
fi

if [ $stage == 3 ];
then
    required="sumner"
    if [ $CLUSTER_NAME == $required ];
    then
        echo "Submitting stage ${stage}"
        sbatch submit_${CLUSTER_NAME}.sb ${stage} ${samplesheet} ${workdir}
    else
        echo "ERORR: To run stage ${stage} ssh to $USER@login.${required}.jax.org"
    fi
fi
