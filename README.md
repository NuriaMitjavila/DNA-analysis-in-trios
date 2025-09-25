# Applied Precision Medicine - DNA analysis in trios
This page is a summary of the report for the project DNA analysis in trios, together with the scripts used to perform the analysis for Case 3: 

"_A medical doctor at medical genetics meets a young girl with frequent coughing. Additionally, she recently suffered a severe lung inflammation. She still has obvious difficulties coughing and breathing, and she is constantly out of breath. Her appearance is that of a child rather small and slim for her age. Neither small stature nor sensitivity to lung inflammation run in the family. The family stems from Southern Europe. The medical doctor suspects a genetic cause and orders whole genome sequencing of the patient and the parents to confirm the diagnosis. The MD asks you to analyze and evaluate the sequencing results._"

For this project, we only focused on *de novo* variants, starting from whole genome sequencing (WGS) data of an affected child and the parents (unaffected). Project data are synthetic with the pathogenic variant mutated manually (find it under the data folder). Here are the steps we followed for the trio analysis: [(0) Getting used to the project directory](#step-0---getting-used-to-the-project-directory), [(1) Read alignment](#step-1---read-alignment), [(2) Variant calling](#step-2---variant-calling), [(3) Joint Genotyping](#step-3---joint-genotyping), [(4) Variant filtering](#step-4---variant-filtering), [(5) De novo detection and more filters](#step-5---de-novo-detection-and-more-filters), [(6) Annotating the de novo candidates](#step-6-annotating-the-de-novo-candidates).<br><br>

### Step 0 - Getting used to the project directory
Download the data from the GitHub repository and store it in the same directory so there is no need to modify the scripts. Here is the distribution: 
- `data` contains the data necessary to perform the analysis: `case3` contains trio WGS data, `ClinVar` contains resources from the database of clinically relevant variants and `reference` contains the human genome assenbly (GRCh38) to ensure all analyses are performed against a standardized coordinate system.
- `scripts` contains the commands that are going to be used to run the analysis with one folder for each step.
- `output` contains all the files that one is creating when running the different scripts. <br><br>

### Step 1 - Read alignment
Raw sequencing reads in FASTQ format were aligned to the human reference genome GRCh38 using the Burrows-Wheeler Aligner. Paired-end reads from each sample were processed in parallel. Intermediate SAM and BAM files were stored in a node-local temporary directory to optimize speed and storage. Aligned reads were converted into CRAM format, indexed, and stored for downstream analysis. **You can find the script in: "scripts/step1"** <br><br>

### Step 2 - Variant calling
Aligned CRAM files were analysed with the Genome Analysis Toolkit (GATK) HaplotypeCaller in GVCF (Genomic Variant Call Format) mode to identify genomic positions differing from the reference. GVCF output captures genotype likelihoods at every position, enabling joint genotyping across multiple samples. HaplotypeCaller was executed using multiple threads to accelerate processing. **You can find the script in: "scripts/step2"** <br><br>

### Step 3 - Joint genotyping
The next step is joint genotyping, which combines the per-sample gVCF into one multi-sample. First, per-sample gVCFs were compressed and indexed using bgzip and tabix to enable efficient random access. Next, multi-sample genotyping was performed with GATK GenomicsDBImport and GenotypeGVCFs, partitioning the genome by chromosome to increase efficiency. Finally, after getting the per-chromosome VCFs, they were merged using GATK MergeVcfs and indexed for further analysis. **You can find the script in: "scripts/step3"** <br><br>

### Step 4 - Variant filtering
Subsequently, variants were filtered to distinguish likely true calls from sequencing artifacts. Raw variant calls often include false positives caused by sequencing errors, misalignments, or technical noise. To address this, we removed low-confidence calls using GATK VariantFiltration with hard thresholds: Quality by Depth (QD < 2.0), Fisher Strand Bias (FS > 60.0), and Mapping Quality (MQ < 40.0). Variant failing any filter were flagged, whereas high-confidence variants were retained for de novo analysis. **You can find the script in: "scripts/step4"** <br><br>

### Step 5 - De novo detection and more filters
The filtered variants were then exported per sample into tabular (TSV) format, including chromosome, position, reference and alternate alleles, genotype (GT), genotype quality (GQ), read depth (DP), and allele depth (AD). In R, the three sample tables were merged on chromosomal position and processed to extract candidate de novo variants, defined as heterozygous in the child and homozygous reference in both patents. Additional quality control included thresholds for GQ (≥20), DP (≥15), allele balance (0.1–0.7), and exclusion of variants with long indels (>10 bp) to focus on high-confidence short variants. **You can find the script in: "scripts/step5"** <br><br>

### Step 6 Annotating the de novo candidates

