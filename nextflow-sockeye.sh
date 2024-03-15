#!/bin/bash

module load git apptainer

if [[ "$USER" == "mbironla" ]]; then
    export CLUSTER_OPTIONS="--account=st-tdjc-1"
else
    export CLUSTER_OPTIONS="--account=st-alexbou-1"
fi

./nextflow $@ -profile cluster
