#!/usr/bin/env bash

#Establishing the CIRCE project directory
CIRCE_DIR="/home/david/BCEMSync/04_Algoritmos/01_proyecto_circRNA/CIRCE";

#Reference genome of the organism
REF_GENOME="/home/david/Homo_sapiens_chr1.fa";

#Input BAM file with RNA-seq data
BAM_FILE="/home/david/chr1_BWA_SRR445016_GRCH38_coord_sorted.bam";

#Log output name
OUTPUT_NAME="Final_verification";

#Running the program
java -d64 -Xmx14g -javaagent:${CIRCE_DIR}/lib/classmexer.jar -cp "${CIRCE_DIR}/lib/NGSEPcore_3.1.1.jar:${CIRCE_DIR}/bin" circe.main.CIRCE ${BAM_FILE} ${REF_GENOME} 1>${CIRCE_DIR}/logs/${OUTPUT_NAME}.circ 2> ${CIRCE_DIR}/logs/${OUTPUT_NAME}.log;