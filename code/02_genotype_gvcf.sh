#!/bin/bash

########### Run as: bash file.sh X 2>/path/to/log.err to run for chrX

echo "Start at $(date)"
set -e

# Tools
GATK=/mnt/data2/saadat/tools/gatk/gatk-4.2.2.0/gatk

# Path
INPUT_DIR=/mnt/data2/saadat/RSV/processed_data/genomics_db_import_output
OUTPUT_DIR=/mnt/data2/saadat/RSV/processed_data/genotype_gvcf_output
REF=/mnt/data2/saadat/Reference_Genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
TEMP_DIR=/mnt/data2/saadat/temp/genotype_gvcf

# Prepare directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

# Prepare input
cd ${INPUT_DIR}

# Run genotype gvcf
$GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms28G -Xmx30G -XX:ParallelGCThreads=2" \
	GenotypeGVCFs \
	-R $REF \
	-V gendb://chr${1}_gdb \
	-O ${OUTPUT_DIR}/chr${1}.g.vcf

echo "End at $(date)"