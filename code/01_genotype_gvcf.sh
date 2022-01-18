#!/bin/bash

########### Run as: bash file.sh 2>/path/to/log.err 1>path/to/log.out

echo "Start at $date"
set -e

# Tools
GATK=/mnt/data2/saadat/tools/gatk/gatk-4.2.2.0/gatk

# Path
INPUT_DIR=/mnt/data2/saadat/RSV/raw_data
OUTPUT_DIR=/mnt/data2/saadat/RSV/processed_data/genotype_gvcf_output
REF=/mnt/data2/saadat/Reference_Genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
TEMP_DIR=/mnt/data2/saadat/temp/genotype_gvcf

# Prepare output
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

# Prepare input
total_input=""
for file in ${INPUT_DIR}/RSV_Antonio/*.vcf.gz
do
    total_input="-V $file ${total_input}"
done
for file in ${INPUT_DIR}/RSV_PRI/*.vcf.gz
do
    total_input="-V $file ${total_input}"
done

# Run genotype gvcf
$GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms28G -Xmx30G -XX:ParallelGCThreads=1" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ${OUTPUT_DIR}/${CHR}_${1}_gdb \
        -R $REF \
        ${total_input} \
        --tmp-dir ${TEMP_DIR} \
        --intervals chr${1} \
        --genomicsdb-shared-posixfs-optimizations true

echo "End at $date"