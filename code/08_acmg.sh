#!/bin/bash 
#SBATCH --chdir /scratch/saadat/RSV
#SBATCH --nodes 1
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 1
#SBATCH --mem 40G
#SBATCH --time 04:00:00
#SBATCH --mail-user=ali.saadat@epfl.ch
#SBATCH -J acmg_RSV
#SBATCH --mail-type=END
#SBATCH -o ./log/acmg_RSV_%J.out # Standard output
#SBATCH -e ./log/acmg_RSV_%J.err # Standard error
echo "START AT $(date)"
set -e

# Path
WORK=/work/gr-fe/saadat/RSV
INPUT_DIR=/work/gr-fe/saadat/RSV/processed_data/pre_annotation_output
OUTPUT_DIR=/work/gr-fe/saadat/RSV/processed_data/tapes_output
TAPES=/work/gr-fe/saadat/tools/tapes/tapes/tapes.py

# Create output folder
mkdir -p ${OUTPUT_DIR}

# Activate python3 virtualn env
source /work/gr-fe/saadat/tools/tapes/python_venv_for_tapes/bin/activate

python3 $TAPES annotate -i ${INPUT_DIR}/bcftools_gatk_vcftools50_norm.vcf.gz  -o ${OUTPUT_DIR}/annotated50.vcf --acmg 
python3 $TAPES sort -i ${OUTPUT_DIR}/annotated50.hg38_multianno.vcf -o ${OUTPUT_DIR}/acmg_results50/ --tab -t 4

echo "END AT $(date)"                                                                                                                                                    1,1           Top