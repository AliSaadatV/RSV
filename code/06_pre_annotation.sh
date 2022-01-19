#!/bin/bash 
#SBATCH --chdir /scratch/saadat/RSV
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 35G
#SBATCH --time 07:00:00
#SBATCH --mail-user=ali.saadat@epfl.ch
#SBATCH -J pre_annotation_RSV
#SBATCH --mail-type=END
#SBATCH -o ./log/pre_annotation_%J.out # Standard output
#SBATCH -e ./log/pre_annotation_%J.err # Standard error
set -e
echo "START AT $(date)"

# Tools
GATK=/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk
BCFTOOLS=/work/gr-fe/saadat/tools/bcftools/installation/bin/bcftools
VCFTOOLS=/work/gr-fe/saadat/tools/vcftools/installation/bin/vcftools
VT=/work/gr-fe/saadat/tools/vt/vt/vt

# Path
WORK=/work/gr-fe/saadat/RSV
SCRATCH=/scratch/saadat/pri/second_try
INPUT_DIR=${WORK}/processed_data/genotype_refinement_output
OUTPUT_DIR=${WORK}/processed_data/pre_annotation_output
REF=/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz

# Create directories
mkdir -p ${SCRATCH}/temp/${SLURM_JOBID}/io
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}
# Remove Batch Effect
$BCFTOOLS filter -i 'QUAL>=25 & INFO/DP>=20 & FORMAT/DP>=10 & FORMAT/GQ>=20' -S . -Ov \
	-o bcftools.vcf --threads 20 \
	${INPUT_DIR}/indel_recal_95_snp_recal_99.7_refined_GQ20.vcf.gz
bgzip -@ 20 -c bcftools.vcf > bcftools.vcf.gz
tabix -p vcf bcftools.vcf.gz
rm bcftools.vcf 

# Exclude filtered
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH}/temp/${SLURM_JOBID}/io -Xms18G -Xmx18G -XX:ParallelGCThreads=4" \
        SelectVariants \
        -V ${OUTPUT_DIR}/bcftools.vcf.gz \
        -R $REF \
        --exclude-filtered \
        -exclude-non-variants \
        --remove-unused-alternates \
        --set-filtered-gt-to-nocall \
        -O ${OUTPUT_DIR}/bcftools_gatk.vcf.gz

# Exclude too variants with minor_AC<=1
$VCFTOOLS --gzvcf bcftools_gatk.vcf.gz --mac 2 --recode --recode-INFO-all --out bcftools_gatk_vcftools
bgzip -@ 20 -c bcftools_gatk_vcftools.recode.vcf > bcftools_gatk_vcftools.vcf.gz
tabix -p vcf bcftools_gatk_vcftools.vcf.gz

# Normalize and break multiallelic 
$VT decompose -s -o temp.vcf bcftools_gatk_vcftools.vcf.gz
$VT normalize -r $REF -o bcftools_gatk_vcftools_norm.vcf temp.vcf
bgzip -@ 20 -c bcftools_gatk_vcftools_norm.vcf > bcftools_gatk_vcftools_norm.vcf.gz
tabix -p vcf bcftools_gatk_vcftools_norm.vcf.gz
rm temp.vcf

echo "FINISH AT $(date)"