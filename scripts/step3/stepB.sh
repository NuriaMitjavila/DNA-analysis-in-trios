#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J joint_genotyping
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-23 

# Load modules
module load bioinfo-tools
module load GATK

# Project paths
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
REF="$PROJECT_FOLDER/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTPUT_FOLDER="$PROJECT_FOLDER/nuria/output"

# Chromosome list (no Y since is a female)
CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHR=${CHRS[${SLURM_ARRAY_TASK_ID}-1]}

# Define workspace and output per chromosome
DB_WORKSPACE="${OUTPUT_FOLDER}/trio_db_${CHR}"
OUT_VCF="${OUTPUT_FOLDER}/trio_joint_chr${CHR}.vcf.gz"

echo "$(date) Chromosome $CHR: Importing gVCFs into GenomicsDB..."

gatk GenomicsDBImport \
    --genomicsdb-workspace-path "$DB_WORKSPACE" \
    -V "${OUTPUT_FOLDER}/child_variants.g.vcf.gz" \
    -V "${OUTPUT_FOLDER}/father_variants.g.vcf.gz" \
    -V "${OUTPUT_FOLDER}/mother_variants.g.vcf.gz" \
    --intervals "$CHR" \
    --reader-threads 16

echo "$(date) Chromosome $CHR: Running joint genotyping..."

gatk GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://${DB_WORKSPACE}" \
    -O "$OUT_VCF"

echo "$(date) Chromosome $CHR done."

