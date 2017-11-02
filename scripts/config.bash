#!/bin/sh -l

# User modifiable variables

PSHOME=/cga/wu/sachet/hla/hla_caller/capture/polysolver_081114
SAMTOOLS_DIR=/cga/wu/sachet/software/samtools
JAVA_DIR=/broad/software/free/Linux/redhat_5_x86_64/pkgs/oracle-java-jdk_1.7.0-17_x86_64/bin
NOVOALIGN_DIR=/broad/software/free/Linux/redhat_5_x86_64/pkgs/novocraft-2.07.18
GATK_DIR=/cga/wu/sachet/hla/hla_caller/capture/polysolver_070914/bin
MUTECT_DIR=/cga/wu/sachet/hla/hla_caller/capture/polysolver_070914/bin
STRELKA_DIR=/cga/wu/sachet/software/bin
TMP_DIR=/cga/wu/sachet
NUM_THREADS=8

export PSHOME
export SCRIPTS_DIR
export SAMTOOLS_DIR
export JAVA_DIR
export NOVOALIGN_DIR
export GATK_DIR
export MUTECT_DIR
export STRELKA_DIR
export TMP_DIR
export NUM_THREADS
