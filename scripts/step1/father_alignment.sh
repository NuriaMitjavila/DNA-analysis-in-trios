#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 16
#SBATCH -t 20:00:00
#SBATCH -J alignment_father
#SBATCH --mail-user=nuriamitjavila2003@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bwa bwa-mem2 samtools

INTERMEDIATE_SAM="$TMPDIR/aligned_reads_father.sam"
INTERMEDIATE_BAM="$TMPDIR/aligned_reads_father.bam"

PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
REF="$PROJECT_FOLDER/reference"
OUTPUT_CRAM="$PROJECT_FOLDER/nuria/output/father_aligned.cram"

bwa-mem2 mem -t 16 -R '@RG\tID:father\tSM:father\tPL:ILLUMINA' \
    $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    $PROJECT_FOLDER/case3/father/forward.fq.gz \
    $PROJECT_FOLDER/case3/father/reverse.fq.gz \
> "$INTERMEDIATE_SAM"

samtools sort -@ 16 -o "$INTERMEDIATE_BAM" "$INTERMEDIATE_SAM"

samtools view -@ 16 -C -T $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -o "$OUTPUT_CRAM" "$INTERMEDIATE_BAM"

samtools index "$OUTPUT_CRAM"

rm "$INTERMEDIATE_SAM" "$INTERMEDIATE_BAM"
