#!/bin/bash 
#SBATCH --chdir /scratch/saadat/RSV
#SBATCH --nodes 1
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 1
#SBATCH --mem 40G
#SBATCH --time 20:00:00
#SBATCH --mail-user=ali.saadat@epfl.ch
#SBATCH -J intervar_RSV
#SBATCH --mail-type=END
#SBATCH -o ./log/intervar_RSV_%J.out # Standard output
#SBATCH -e ./log/intervar_RSV_%J.err # Standard error
echo "START AT $(date)"
set -e

# Path
WORK=/work/gr-fe/saadat/RSV
INPUT_DIR=/work/gr-fe/saadat/RSV/processed_data/pre_annotation_output
OUTPUT_DIR=/work/gr-fe/saadat/RSV/processed_data/intervar_output
INTERVAR=/work/gr-fe/saadat/tools/inter_var/InterVar-master
ANNOVAR=/work/gr-fe/saadat/tools/annovar/annovar

# Create output folder
mkdir -p ${OUTPUT_DIR}

# Activate load python3
module load gcc python/3.7.7

python ${INTERVAR}/Intervar.py -i ${INPUT_DIR}/bcftools_gatk_vcftools50_norm.vcf.gz \
    -o ${OUTPUT_DIR}/intervar50.vcf \
    -b hg38 \
    --input_type=VCF_m \
    -t ${INTERVAR}/intervardb \
    --table_annovar=${INTERVAR}/table_annovar.pl \
    --convert2annovar=${INTERVAR}/convert2annovar.pl \
    --annotate_variation=${INTERVAR}/annotate_variation.pl \
    -d ${ANNOVAR}/humandb_hg38:

echo "END AT $(date)"                                                                                                                                                    1,1           Top