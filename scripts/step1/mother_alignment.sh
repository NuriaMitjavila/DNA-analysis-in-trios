#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 16
#SBATCH -t 20:00:00
#SBATCH -J alignment_mother
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bwa bwa-mem2 samtools

INTERMEDIATE_SAM="$TMPDIR/aligned_reads_mother.sam"
INTERMEDIATE_BAM="$TMPDIR/aligned_reads_mother.bam"

PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/private/nurmi143/project"
REF="$PROJECT_FOLDER/data/reference"
OUTPUT_CRAM="$PROJECT_FOLDER/output/mother_aligned.cram"

bwa-mem2 mem -t 16 -R '@RG\tID:mother\tSM:mother\tPL:ILLUMINA' \
    $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    $PROJECT_FOLDER/data/case3/mother/forward.fq.gz \
    $PROJECT_FOLDER/data/case3/mother/reverse.fq.gz \
> "$INTERMEDIATE_SAM"

samtools sort -@ 16 -o "$INTERMEDIATE_BAM" "$INTERMEDIATE_SAM"

samtools view -@ 16 -C -T $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -o "$OUTPUT_CRAM" "$INTERMEDIATE_BAM"

samtools index "$OUTPUT_CRAM"

rm "$INTERMEDIATE_SAM" "$INTERMEDIATE_BAM"
