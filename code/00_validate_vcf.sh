#!/bin/bash

########### Run as: bash file.sh 2>>/path/to/log.err 1>>path/to/log.out

echo "Start at $date"
set -e

# Tools
GATK=/mnt/data2/saadat/tools/gatk/gatk-4.2.2.0/gatk

# Path
INPUT_DIR=/mnt/data2/saadat/RSV/raw_data/RSV_Antonio
REF=/mnt/data2/saadat/Reference_Genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz

# Validate VCF
for file in ${INPUT_DIR}/*.vcf.gz
do
    $GATK ValidateVariants \
    -R $REF \
    -V $file \
    --validation-type-to-exclude ALL
done

echo "End at $date"