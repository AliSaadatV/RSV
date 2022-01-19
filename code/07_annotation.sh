#!/bin/bash 
#SBATCH --chdir /scratch/saadat/RSV
#SBATCH --nodes 1
#SBATCH --ntasks 10
#SBATCH --cpus-per-task 1
#SBATCH --mem 30G
#SBATCH --time 03:00:00
#SBATCH --mail-user=ali.saadat@epfl.ch
#SBATCH -J anno_RSV
#SBATCH --mail-type=END
#SBATCH -o ./log/anno_RSV_%J.out # Standard output
#SBATCH -e ./log/anno_RSV_%J.err # Standard error
echo "START AT $(date)"
set -e

# Path
WORK=/work/gr-fe/saadat/RSV
REF=/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
SCRATCH=/scratch/saadat/RSV
CACHE=/work/gr-fe/databases/vep_hg38/cache_hg38_ensembl
PLUGINS=/work/gr-fe/databases/vep_hg38/plugins
INPUT_DIR=/work/gr-fe/saadat/RSV/pre_annotation_output
OUTPUT_DIR=/work/gr-fe/saadat/RSV/annotation_output
LOFTEE=/work/gr-fe/databases/vep_hg38/loftee
ANNOVAR=/work/gr-fe/saadat/tools/annovar/annovar
SLIVAR_DIR=/work/gr-fe/saadat/tools/slivar

# Create output folder
mkdir -p ${OUTPUT_DIR}

# Activate vep from conda
eval "$(conda shell.bash hook)"
conda activate ensembl-vep

vep \
--offline --vcf \
--fasta $REF \
--cache \
--dir_cache $CACHE \
--canonical \
--hgvs \
--symbol \
--mane \
--flag_pick \
--cache_version 104 \
--dir_plugins ${LOFTEE}/plugin_hg38/loftee \
--plugin LoF,loftee_path:${LOFTEE}/plugin_hg38/loftee,human_ancestor_fa:${LOFTEE}/data_hg38/human_ancestor.fa.gz,gerp_bigwig:${LOFTEE}/data_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
-i ${INPUT_DIR}/bcftools_gatk_vcftools_norm.vcf.gz \
-o ${OUTPUT_DIR}/vep.vcf.gz \
--stats_file ${OUTPUT_DIR}/vep.html \
--fork 10 \
--force_overwrite

# Filter with Slivar
${SLIVAR_DIR}/slivar expr \
--js ${SLIVAR_DIR}/slivar-functions.js \
--vcf ${OUTPUT_DIR}/vep.vcf.gz \
--info 'INFO.impactful && variant.FILTER == "PASS"' \
--out-vcf ${OUTPUT_DIR}/vep_slivar.vcf

# Annotate with annovar
perl ${ANNOVAR}/table_annovar.pl ${OUTPUT_DIR}/vep_slivar.vcf ${ANNOVAR}/humandb_hg38/ \
 -buildver hg38 \
 -thread 10 \
 -out ${OUTPUT_DIR}/vep_slivar_annovar.vcf \
 -vcfinput \
 -remove \
 -protocol refGene,ensGene,dbnsfp42a,dbnsfp31a_interpro,avsnp147,clinvar_20210501,intervar_20180118,esp6500siv2_all,gnomad30_genome,exac03,ALL.sites.2015_08,dbscsnv11 \
 -operation g,g,f,f,f,f,f,f,f,f,f,f \
 -nastring .

echo "END AT $(date)"                                                                                                                                                    1,1           Top