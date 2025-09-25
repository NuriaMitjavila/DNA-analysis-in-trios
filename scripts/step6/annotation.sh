# gets the first 4 columns: CHROM, POS, REF, ALT
cut -f1-4 de_novo_candidates_filtered.tsv > de_novo_for_annotation.tsv

# create a BED file with variant positions: VCF uses 1-based coordinates, but BED files require 0-based start and 1-based end coordinates, so $2-1 gives the start, $2 gives the end.
awk '{print $1"\t"$2-1"\t"$2}' de_novo_for_annotation.tsv > positions.bed

# Remove the header line: +2 means “start from line 2 onward”
tail -n +2 positions.bed > positions_noheader.bed

# Extract variants from the VCF using positions
FILTERED_VCF="$PROJECT_FOLDER/nuria/output/trio_joint.filtered.vcf.gz"
bcftools view -R positions_noheader.bed $FILTERED_VCF -o de_novo_candidates.vcf -O z

# Compress and index the candidate VCF
bgzip de_novo_candidates.vcf
tabix -p vcf de_novo_candidates.vcf.gz
