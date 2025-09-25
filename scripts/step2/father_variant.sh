#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J variantcalling_father
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL

# Load necessary modules
module load bioinfo-tools
module load GATK

# Project paths
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
REF="$PROJECT_FOLDER/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT_CRAM="$PROJECT_FOLDER/nuria/output/father_aligned.cram"
OUTPUT_GVCF="$PROJECT_FOLDER/nuria/output/father_variants.g.vcf"

# Run GATK HaplotypeCaller in GVCF mode
gatk HaplotypeCaller \
    -R $REF \
    -I $INPUT_CRAM \
    -O $OUTPUT_GVCF \
    -ERC GVCF \
    --native-pair-hmm-threads 8

