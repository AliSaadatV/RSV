#!/bin/bash

########### Run as: bash file.sh 

echo "Start at $(date)"
set -e

# Tools
GATK=/mnt/data2/saadat/tools/gatk/gatk-4.2.2.0/gatk

# Path
INPUT_DIR=/mnt/data2/saadat/RSV/processed_data/genotype_gvcf_output
OUTPUT_DIR=/mnt/data2/saadat/RSV/processed_data/gather_vcfs_output
TEMP_DIR=/mnt/data2/saadat/temp/gather_vcfs

# Prepare directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

# Prepare input
cd ${INPUT_DIR}
input_list=$( for i in {1..22} X Y M; do echo -n "-I chr${i}.g.vcf " ; done)

# Run genotype gvcf
$GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms28G -Xmx30G -XX:ParallelGCThreads=2" \
	GatherVcfs \
	${input_list}\
	-O ${OUTPUT_DIR}/merged.vcf.gz

echo "End at $(date)"