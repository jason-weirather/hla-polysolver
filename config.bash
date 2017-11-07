#!/bin/sh -l

# User modifiable variables

PSHOME=/home/hla-polysolver
GATK_DIR=$CONDA_PREFIX/jar
MUTECT_DIR=$CONDA_PREFIX/jar
STRELKA_DIR=$CONDA_PREFIX/bin
TMP_DIR=/tmp
NUM_THREADS=8

export PSHOME
export GATK_DIR
export MUTECT_DIR
export STRELKA_DIR
export TMP_DIR
export NUM_THREADS
