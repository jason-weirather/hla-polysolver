#!/bin/sh -l

# User modifiable variables

#export PSHOME=$CONDA_PREFIX/scripts
export GATK_DIR=$CONDA_PREFIX/jar
export MUTECT_DIR=$CONDA_PREFIX/jar
export STRELKA_DIR=$CONDA_PREFIX/include/strelka-upstream-v1.0.11
export DATA_DIR=$CONDA_PREFIX/share/polysolver_data
export TMP_DIR=/tmp
export NUM_THREADS=8
