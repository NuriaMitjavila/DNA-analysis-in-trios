# Applied Precision Medicine - DNA analysis in trios
This page is a summary of the report for the project DNA analysis in trios, together with the scripts used to perform the analysis for Case 3: 

"_A medical doctor at medical genetics meets a young girl with frequent coughing. Additionally, she recently suffered a severe lung inflammation. She still has obvious difficulties coughing and breathing, and she is constantly out of breath. Her appearance is that of a child rather small and slim for her age. Neither small stature nor sensitivity to lung inflammation run in the family. The family stems from Southern Europe. The medical doctor suspects a genetic cause and orders whole genome sequencing of the patient and the parents to confirm the diagnosis. The MD asks you to analyze and evaluate the sequencing results._"

For this project, we only focused on *de novo* variants, starting from whole genome sequencing (WGS) data of an affected child and the parents (unaffected). Project data are synthetic with the pathogenic variant mutated manually (find it under the data folder). Here are the steps we followed for the trio analysis: [(0) Getting used to the project directory](#step-0---getting-used-to-the-project-directory), [(1) Read alignment](#step-1---read-alignment), [(2) Variant calling](#step-2---variant-calling), [(3) Joint Genotyping](#step-3---joint-genotyping), [(4) Variant filtering](#step-4---variant-filtering), [(5) De novo detection and more filters](#step-5---de-novo-detection-and-more-filters), [(6) Annotating the de novo candidates](#step-6-annotating-the-de-novo-candidates).<br><br>

### Step 0 - Getting used to the project directory
From the terminal, use `ssh username@rackham.uppmax.uu.se` and enter your UPPMAX password, you are now in your login node `/home/username`. All project data are stored in the directory `/crex/proj/uppmax2024-2-1/rare_variants` in the following way: 
- `case3` contains trio WGS data with each sample's raw sequencing from high-throughput sequencing platforms (e.g., Illumina).
- `ClinVar` contains resources from the database of clinically relevant variants to provide information on pathogenicity and associated diseases.
- `reference` contains the human genome assenbly (GRCh38) for read alignment and variant calling, to ensure all analyses are performed against a standardized coordinate system.
- `PROJECT_FOLDER` will be the folder where we will be performing the calculations and where the scripts will be allocated.<br><br>

### Step 1 - Read alignment
Raw sequencing reads in FASTQ format were aligned to the human reference genome GRCh38 using the Burrows-Wheeler Aligner. Paired-end reads from each sample were processed in parallel. Intermediate SAM and BAM files were stored in a node-local temporary directory to optimize speed and storage. Aligned reads were converted into CRAM format, indexed, and stored for downstream analysis. **You can find the script in: "scripts/step1"**

## Step 2 - Variant calling
Aligned CRAM files were analysed with the Genome Analysis Toolkit (GATK) HaplotypeCaller in GVCF (Genomic Variant Call Format) mode to identify genomic positions differing from the reference. GVCF output captures genotype likelihoods at every position, enabling joint genotyping across multiple samples. HaplotypeCaller was executed using multiple threads to accelerate processing. **You can find the script in: "scripts/step2"**

## Step 3 - Joint genotyping
The next step is joint genotyping, which combines the per-sample gVCF into one multi-sample. First, per-sample gVCFs were compressed and indexed using bgzip and tabix to enable efficient random access. Next, multi-sample genotyping was performed with GATK GenomicsDBImport and GenotypeGVCFs, partitioning the genome by chromosome to increase efficiency. Finally, after getting the per-chromosome VCFs, they were merged using GATK MergeVcfs and indexed for further analysis. **You can find the script in: "scripts/step3"**

## Step 4 - Variant filtering
Subsequently, variants were filtered to distinguish likely true calls from sequencing artifacts. Raw variant calls often include false positives caused by sequencing errors, misalignments, or technical noise. To address this, we removed low-confidence calls using GATK VariantFiltration with hard thresholds: Quality by Depth (QD < 2.0), Fisher Strand Bias (FS > 60.0), and Mapping Quality (MQ < 40.0). Variant failing any filter were flagged, whereas high-confidence variants were retained for de novo analysis. **You can find the script in: "scripts/step4"**

## Step 5 - De novo detection and more filters

If it isn't specified that we are looking for are *de novo* variants, we would do annotation first. However, annotating all filtered variants takes much longer time and is likely unnecesary. Although tools like `bcftools` can also filter based on the genotypes, we decide to use R for a more visually pleasing and convenient experience.

First, let's quickly check how the genotypes have been encoded in the VCF file, using `bcftools`:

```
module load bioinfo-tools bcftools
bcftools query -f '%CHROM\t%POS\t[%SAMPLE=%GT\t]\n' $FILTERED_VCF | head
```

Here `bcftools query` extracts information from a VCF file in a tabular format according to the format string `-f`, followed by which we define what fields and how we want to print out:
- `%CHROM`: the chromosome of the variant
- `%POS`: the position on the chromosome
- `[%SAMPLE=%GT\t]`: a per-sample loop. For every sample at this site: `%SAMPLE` = the sample name (e.g., child, father, mother), `%GT` = the genotype; `\t` will add a tab in between
- `\n`, just like `\t`, is an escape sequence and stands for the newline character

With the command above, in an interactive window, you can see something like:

```
1	10125	child=./. father=0/1  mother=0/0
1	10439	child=0/0	father=0|1	mother=0/0
1	10583	child=0/0	father=0/0	mother=0/1
1	13613	child=1/1	father=0/1	mother=0/1
1	13649	child=0/0	father=0/1	mother=0/0
1	13684	child=./.	father=0/1	mother=0/1
1	13757	child=0/0	father=0/1	mother=0/0
1	13813	child=0|1	father=0|1	mother=0|1
1	13838	child=0|1	father=0|1	mother=0|1
1	13868	child=0/0	father=0/1	mother=0/1
```

The genotypes are written with numbers separated by `/` or `|` where: 0 = the reference allele (from the reference genome), 1 = the first alternate allele listed in the VCF. So the codes mean: `0/0` = homozygous reference, `0/1` heterozygous, `1/1` homozygous alternate; ./. = missing genotype (no call made for this sample). In the trio example, you would expect the child as `0/1` (`0|1`, if phased) while both parents are `0/0` (`0|0`) for a de novo variant.

For the ease of usage in R, we can export the VCF in to TSV file for each sample (just so we don't mix the field names). Remember to modify the sample fields in the following command:

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%GQ\t%DP\t%AD]\n' --samples sample1 $FILTERED_VCF > $FILTERED_VCF.sample1.tsv
```

- `REF` = reference allele; `ALT` = alternative allele
- `FILTER`, the field generated by `VariantFiltration`
- `GQ` = Genotype quality - higher values mean more confidence that the reported genotype is correct (e.g., GQ=99 means very high confidence)
- `DP` = Read depth - the total number of sequencing reads covering that position for the sample (e.g. DP=35 means 35 reads overlapped this variant site)
- `AD` = Allele depths - a comma-separated count of reads supporting each allele (e.g. "AD=28,7" means 28 reads support the reference allele and 7 reads support the alternate allele); useful for checking allele balance (whether the variant allele is present in a reasonable fraction of reads)

Together these fields can help filter out false positives. For instance, a genotype being 0/1 but AD=30,1 (only 1 alt read) and low GQ might be an artifact, but a de novo candidate with 0/1, AD=15,12, DP=27, and high GQ looks much more convincing. Before continuing, you can extract those that passed the hard filtering using `awk`:

```
awk '$5=="PASS"' $FILTERED_VCF.sample1.tsv > $FILTERED_VCF.sample1.tsv
```

Then, load the R packages by `module load R_packages` and enter `R` to start the R session. Some examples and hints are shown below:

#### Reading and processing the genotype TSV files

For manipulating the dataframe, it is recommended to load the R packages `dplyr` and `tidyr` first. For each sample, you can edit the column names to avoid ambiguity.

```
sample1 <- read.table("sample1.tsv", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(sample1) <- c("CHROM","POS","REF","ALT","sample1_GT","sample1_GQ","sample1_DP","sample1_AD")
```

Since we want to merge the three samples, we need a unique "key" for them to match for. In this case, the chromosome positions:

```
sample1 <- sample1 %>% mutate(KEY = paste(CHROM, POS, sep = ":"))
```

And merge by `inner_join()` (add new lines for the 3rd sample):

```
trio <- sample1 %>%
  inner_join(sample2 %>% select(KEY, sample2_GT = sample2_GT, sample2_AD = sample2_AD, sample2_DP = sample2_DP, sample2_GQ = sample2_GQ), by = "KEY")
```

Remember for allele depths, in one entry they contain two values seperated by a comma, so splitting them into ref/alt will be helpful. For example:

```
trio <- trio %>%
  separate(sample1_AD, into = c("sample1_AD_ref","sample1_AD_alt"), sep = ",", convert = TRUE) %>%
  separate(sample2_AD, into = c("sample2_AD_ref","sample2_AD_alt"), sep = ",", convert = TRUE)
```

You should also check the data types of the fields that should be numeric and convert them (here `ends_with()` looks for the columns with names ending with these strings):

```
trio <- trio %>%
  mutate(across(ends_with(c("AD_ref","AD_alt","DP","GQ")), as.numeric))
```

#### Detecting de novo

Now you can filter for the de novo candidates:

```
de_novo <- trio %>%
  filter(
    child_GT %in% c("0/1","1/0","0|1","1|0"),
    ...... # what should you write for the father and mother's GTs?
  )
```

#### Extra filters

Some additional filters to consider about (adjust thresholds if needed):
- GQ ≥ 20~30
- DP ≥15–20
- Parents' AD_alt == 0  (or ≤1 read with AD_alt very low) plus good DP and GQ
- Allele balance: add a column to the dataframe with `mutate()` as above, and filter for the value among a range, such as 0.3~0.7. Hint: allele balance = AD_alt / (AD_ref + AD_alt)
- For typical trio-based de novo variant discovery, the "true" de novo variants are usually short variants such as SNPs (single-base changes) or small indels (usually ≤10 bp). In the genotype file, some have very long alleles which can usually be structural variants, complex events, or annotation artifacts. Most trio studies would focus on short, high-confidence variants because the long ones are harder to genotype accurately and are rarely de novo. Hint: filter the length of the alleles (ref and alt) using `nchar()`.

Run the filtering criteria one by one, and then output them to a new TSV, e.g.:

```
write.table(de_novo_filtered, "de_novo_candidates_filtered.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

## Step 6 Annotating the de novo candidates

Raw variants by themselves are just genomic coordinates and alleles. Annotation enriches them with biological meaning: predicting their effect on genes (via tools like SnpEff or VEP), linking them to known databases (e.g., dbSNP for rsIDs, ClinVar for disease associations). This step transforms a simple variant list into interpretable results with clinical and research relevance.

We will use SnpEff and ClinVar. SnpEff provides detailed functional predictions, such as whether a variant causes a missense or nonsense change, while ClinVar offers clinically curated information about known pathogenic and benign variants. In our workflow, we did not include dbSNP separately, because ClinVar already incorporates the dbSNP reference identifiers (rsIDs) alongside its clinical interpretations. This makes ClinVar sufficient for both functional and clinical annotation in this teaching context, while keeping the pipeline simpler.

Since we have extracted the de novo candidates into TSV file for inspection and filtering, we need to convert them back into a VCF file before annotation:

```
# gets the first 4 columns: CHROM, POS, REF, ALT
cut -f1-4 de_novo_candidates_filtered.tsv > de_novo_for_annotation.tsv
# create a BED file with variant positions: VCF uses 1-based coordinates, but BED files require 0-based start and 1-based end coordinates, so $2-1 gives the start, $2 gives the end.
awk '{print $1"\t"$2-1"\t"$2}' de_novo_for_annotation.tsv > positions.bed
# Remove the header line: +2 means “start from line 2 onward”
tail -n +2 positions.bed > positions_noheader.bed
# Extract variants from the VCF using positions
bcftools view -R positions_noheader.bed $FILTERED_VCF -o de_novo_candidates.vcf -O z
# Compress and index the candidate VCF
bgzip de_novo_candidates.vcf
tabix -p vcf de_novo_candidates.vcf.gz
```

#### Using SnpEff

Load the module `snpEff` once you have loaded `bioinfo-tools`. SnpEff is written in Java, so it requires a Java runtime environment to execute. On HPC systems, different versions of Java may be available, and not all software is guaranteed to work with every version. We can check the versions of Java available on Rackham by `module avail java`; by `module load java/OpenJDK_17+35`, we explicitly load Java version 17 (OpenJDK build 35), which is known to be compatible with SnpEff.

Then SnpEff shall be able to run smoothly and you can also get the compressed output together:

```
snpEff -Xmx24g -v GRCh38.105 de_novo_candidates.vcf \
  | bcftools view -Oz -o de_novo_snpEff.vcf.gz

tabix -p vcf de_novo_snpEff.vcf.gz
```

- `-Xmx24g`: tells Java to allow SnpEff to use up to 24 GB of memory. Adjust if you need to.
- `-v`: enables verbose mode (prints more details about what SnpEff is doing). `GRCh38.105`: specifies the genome database SnpEff should use for annotation (here GRCh38, Ensembl release 105).

Similarly, for simplicity, you can extract the information and convert them into TSV file:

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' de_novo_snpEff.vcf.gz > de_novo_snpEff.tsv
```

You can either load it in R on Rackham or download the file to your local computer. In short, SnpEff can predict whether the change occurs in a coding region, whether it alters the protein sequence, and what gene/transcript is affected. These predictions are stored in the `ANN` field of the `INFO` column in the annotated VCF, which you can access by `INFO/ANN`. For more details, you can check by:

```
bcftools view -h de_novo_snpEff.vcf.gz | grep "^##INFO"
```

We are interested in the impact categories:
- HIGH: Likely to have a strong effect on the protein (e.g., stop-gain, frameshift, splice donor/acceptor changes).
- MODERATE: May change protein function but not as drastically (e.g., missense variants that swap one amino acid for another).
- LOW: Less likely to alter protein function (e.g., synonymous changes that don’t change the amino acid).
- MODIFIER: Variants outside coding regions or with unknown impact (e.g., intergenic or intronic variants).

To filter for candidates that are more biologically meaningful and likely to contribute to disease, we can focus on "HIGH" and "MODERATE".

#### Using ClinVar

ClinVar is a large public database that links genetic variants to clinical interpretations, such as whether a variant is benign, likely benign, pathogenic, or of uncertain significance. When we annotate our candidate variants with ClinVar and then inspect them using `bcftools`. For example:

```
bcftools annotate -a /crex/proj/uppmax2024-2-1/rare_variants/ClinVar/clinvar_20250831.vcf.gz -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/ORIGIN,INFO/RS \
    de_novo_candidates.vcf.gz -Oz -o de_novo_ClinVar.vcf.gz
tabix -p vcf de_novo_ClinVar.vcf.gz
```

Before exporting your results, you may check what annotation has ClinVar given by the `bcftools view` command shown above. Typically, we would be interested in the clinical significance (e.g., "pathogenic"), how strong the supporting evidence is (e.g., criteria provided, multiple submitted, no conflicts), the disease or condition associated with the variant, any database identifiers (e.g., OMIM, MedGen) linked to that condition, the dbSNP identifier if available. In addition, we may be curious about the allele origin, which is the type of sample the variant was observed in. Instead of storing plain text like germline or somatic, ClinVar compresses them into a bit flag integer. Each bit in the number corresponds to one possible origin. If you check ClinVar's documentaion or the header information of your annotated output, you may see:

```
##INFO=<ID=ORIGIN,Number=1,Type=Integer,
Description="Allele origin. Bit-encoded: 
1-germline, 2-somatic, 4-inherited, 8-paternal, 
16-maternal, 32-de-novo, 64-biparental, 
128-uniparental, 256-not-tested, 512-tested-inconclusive, 
1073741824-other">
```

Now you know what to look for and extract!
