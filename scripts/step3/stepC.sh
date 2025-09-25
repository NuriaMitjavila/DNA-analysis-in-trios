#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J merge_vcfs
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools
module load GATK

# Project paths
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
OUTPUT_FOLDER="$PROJECT_FOLDER/nuria/output"

# Gather all per-chromosome VCFs
INPUTS=(${OUTPUT_FOLDER}/trio_joint_chr*.vcf.gz)

# Merge into one genome-wide VCF
gatk MergeVcfs $(printf -- "-I %s " "${INPUTS[@]}") -O "${OUTPUT_FOLDER}/trio_joint.vcf.gz"

# Index the merged VCF
gatk IndexFeatureFile -I "${OUTPUT_FOLDER}/trio_joint.vcf.gz"

echo "Genome-wide joint VCF created: ${OUTPUT_FOLDER}/trio_joint.vcf.gz"

