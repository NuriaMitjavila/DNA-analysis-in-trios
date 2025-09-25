module load bioinfo-tools bcftools

FILTERED_VCF="/crex/proj/uppmax2024-2-1/rare_variants/nuria/output/trio_joint.filtered.vcf.gz"

# Extract info per sample
for SAMPLE in child father mother; do
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%GQ\t%DP\t%AD]\n' \
        --samples $SAMPLE \
        $FILTERED_VCF > "${FILTERED_VCF}.${SAMPLE}.tsv"

    # Keep only variants that PASS
    awk '$5=="PASS"' "${FILTERED_VCF}.${SAMPLE}.tsv" > "${FILTERED_VCF}.${SAMPLE}.tsv.pass"
done

