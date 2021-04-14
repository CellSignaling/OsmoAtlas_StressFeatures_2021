###################
#
# 
# Script for gathering genes from HeLa analysis
# 
#
###################
library(dplyr)

## Experimental
print("Gathering experimental data")
# HeLa
resHeLa = read.delim("data/original_data/Cate_DESeq2_salt_vs_nostress.tsv", stringsAsFactors = F)

# Remove those genes with padj(salt comparison) == NA, because we can't conclude that they are unresponsive
resHeLa = resHeLa[!is.na(resHeLa$padj),]
resHeLa = resHeLa %>%
  mutate(salt_HeLa = ifelse(log2FoldChange >= 1 & padj <= 0.05, yes = "upregulated", no = ifelse(log2FoldChange <= -1 & padj <= 0.05, yes = "downregulated", no = "unresponsive")))


# Format output
resTable.o = resHeLa %>%
  select("X", "gene", "salt_HeLa")

colnames(resTable.o)[1:2] = c("ensembl_gene", "associated_gene_name")

write.table(x = resTable.o,  "results/human/humanStressConsensusTable_experimental_HeLa.tab", quote = F, row.names = F, sep = "\t")
