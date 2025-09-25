library(dplyr)
library(tidyr)

# Read per-sample TSVs
setwd("C:/Users/nuria/Downloads/")
child <- read.table("child_results.tsv", header=FALSE, sep="\t", stringsAsFactors=FALSE)
father <- read.table("father_results.tsv", header=FALSE, sep="\t", stringsAsFactors=FALSE)
mother <- read.table("mother_results.tsv", header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Rename columns
colnames(child) <- c("CHROM","POS","REF","ALT","FILTER","child_GT","child_GQ","child_DP","child_AD")
colnames(father) <- c("CHROM","POS","REF","ALT","FILTER","father_GT","father_GQ","father_DP","father_AD")
colnames(mother) <- c("CHROM","POS","REF","ALT","FILTER","mother_GT","mother_GQ","mother_DP","mother_AD")

# Add KEY for merging
child <- child %>% mutate(KEY = paste(CHROM, POS, sep=":"))
father <- father %>% mutate(KEY = paste(CHROM, POS, sep=":"))
mother <- mother %>% mutate(KEY = paste(CHROM, POS, sep=":"))

# Merge trio
trio <- child %>%
  inner_join(father %>% select(KEY, father_GT, father_GQ, father_DP, father_AD), by="KEY") %>%
  inner_join(mother %>% select(KEY, mother_GT, mother_GQ, mother_DP, mother_AD), by="KEY")

# Split ref/alt counts:
trio <- trio %>%
  separate(child_AD, into=c("child_AD_ref","child_AD_alt"), sep=",", convert=TRUE) %>%
  separate(father_AD, into=c("father_AD_ref","father_AD_alt"), sep=",", convert=TRUE) %>%
  separate(mother_AD, into=c("mother_AD_ref","mother_AD_alt"), sep=",", convert=TRUE)

# Check the data types that should be numeric
trio <- trio %>% mutate(across(ends_with(c("AD_ref","AD_alt","DP","GQ")), as.numeric))

# Detecting de novo
de_novo <- trio %>% filter(child_GT %in% c("0/1","1/0","0|1","1|0"), father_GT %in% c("0/0","0|0"), mother_GT %in% c("0/0","0|0"))

# Extra quality filters
de_novo_filtered <- de_novo %>%
   mutate(child_AB = child_AD_alt / (child_AD_ref + child_AD_alt), ref_len = nchar(REF), alt_len = nchar(ALT)) %>%
   filter(child_GQ >= 30, child_DP >= 20, father_AD_alt == 0, mother_AD_alt == 0, child_AB >= 0.1 & child_AB <= 0.7, ref_len <= 10, alt_len <= 10)

# Save the results
write.table(de_novo_filtered, "de_novo_candidates_filtered.tsv", sep="\t", quote=FALSE, row.names=FALSE)

