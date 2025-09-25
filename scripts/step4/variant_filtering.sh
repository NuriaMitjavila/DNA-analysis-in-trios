#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J variant_filtering
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools
module load GATK
module load bcftools

# Project paths
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
REF="$PROJECT_FOLDER/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT_VCF="$PROJECT_FOLDER/nuria/output/trio_joint.vcf.gz"
FILTERED_VCF="$PROJECT_FOLDER/nuria/output/trio_joint.filtered.vcf.gz"

echo "$(date) Applying hard filters to joint VCF..."

gatk VariantFiltration \
   -R $REF \
   -V $INPUT_VCF \
   -O $FILTERED_VCF \
   --filter-name "QD_filter" --filter-expression "QD < 2.0" \
   --filter-name "FS_filter" --filter-expression "FS > 60.0" \
   --filter-name "MQ_filter" --filter-expression "MQ < 40.0"

echo "$(date) Indexing filtered VCF..."

bcftools index -t $FILTERED_VCF

echo "Filtering complete. Filtered VCF: $FILTERED_VCF"

