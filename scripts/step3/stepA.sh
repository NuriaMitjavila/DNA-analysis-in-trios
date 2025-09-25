#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J compress_gvcfs
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL

# Load modules
module load bioinfo-tools
module load htslib

# Project paths
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
OUTPUT_FOLDER="$PROJECT_FOLDER/nuria/output"

# List of samples
SAMPLES=("child" "father" "mother")

# Loop through samples
for SAMPLE in "${SAMPLES[@]}"; do
    INPUT_GVCF="$OUTPUT_FOLDER/${SAMPLE}_variants.g.vcf"
    COMPRESSED_GVCF="$OUTPUT_FOLDER/${SAMPLE}_variants.g.vcf.gz"

    if [ ! -f "$COMPRESSED_GVCF" ]; then
        echo "Compressing $SAMPLE..."
        bgzip -c "$INPUT_GVCF" > "$COMPRESSED_GVCF"
        tabix -p vcf "$COMPRESSED_GVCF"
    else
        echo "$SAMPLE gVCF already compressed and indexed."
    fi
done

echo "All gVCFs compressed and indexed."

