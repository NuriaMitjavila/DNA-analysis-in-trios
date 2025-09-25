bcftools annotate -a /crex/proj/uppmax2024-2-1/rare_variants/ClinVar/clinvar_20250831.vcf.gz -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/ORIGIN,INFO/RS de_novo_candidates.vcf.gz -Oz -o de_novo_ClinVar.vcf.gz

tabix -p vcf de_novo_ClinVar.vcf.gz
