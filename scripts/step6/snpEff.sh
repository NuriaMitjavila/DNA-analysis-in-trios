module load snpEff
module load java/OpenJDK_17+35

snpEff -Xmx24g -v GRCh38.105 de_novo_candidates.vcf | bcftools view -Oz -o de_novo_snpEff.vcf.gz
tabix -p vcf de_novo_snpEff.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' de_novo_snpEff.vcf.gz > de_novo_snpEff.tsv
cat de_novo_snpEff.tsv | grep "HIGH" > de_novo_snpEff_impact.tsv
cat de_novo_snpEff.tsv | grep "MODERATE" >> de_novo_snpEff_impact.tsv
