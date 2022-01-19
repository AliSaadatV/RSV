#!/bin/bash

########### Run as: bash file.sh 

echo "Start at $(date)"
set -e

# Tools
GATK=/mnt/data2/saadat/tools/gatk/gatk-4.2.2.0/gatk

# Path
INPUT_DIR=/mnt/data2/saadat/RSV/processed_data/vqsr_output
OUTPUT_DIR=/mnt/data2/saadat/RSV/processed_data/genotype_refinement_output
REF=/mnt/data2/saadat/Reference_Genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
TEMP_DIR=/mnt/data2/saadat/temp/genotype_refinement
KNOWN_SITES=/mnt/data2/saadat/datasets

# Prepare directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

# Run genotype refinement
$GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms25G -Xmx30G -XX:ParallelGCThreads=2" \
	CalculateGenotypePosteriors \
	-V ${INPUT_DIR}/indel_recal_95_snp_recal_99.7.vcf.gz \
	-O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined.vcf.gz \
	--supporting-callsets ${KNOWN_SITES}/af-only-gnomad.hg38.vcf.gz \
	--num-reference-samples-if-no-call 20314 && \

$GATK --java-options "-Djava.io.tmpdir=${TEMP_DIR} -Xms25G -Xmx30G -XX:ParallelGCThreads=2" \
    VariantFiltration \
    -V ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined.vcf.gz \
    -R $REF \
    --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "lowGQ" \
    -O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined_GQ20.vcf.gz

echo "End at $(date)"