#!/bin/sh -l

# User modifiable variables

export PSHOME=/home/hla-polysolver
export GATK_DIR=$CONDA_PREFIX/jar
export MUTECT_DIR=$CONDA_PREFIX/jar
export STRELKA_DIR=$CONDA_PREFIX/bin
export DATA_DIR=$CONDA_PREFIX/share/polysolver_data
export TMP_DIR=/tmp
export NUM_THREADS=8
