#################################
#
# 
# TCGA CNA analysis
# 
# 
#################################

## Load required libraries
library(pbapply)
library(plyr)

tumours = dir("data/original_data/tcga.cn", full.names = T)


cnProp.lst = pblapply(tumours, cl = 4, function(tumour){
  tumourName = basename(tumour)
  print(paste("Computing CNA proportions from the following tumor type: ----->", tumourName))
  tumourFile = paste0(tumour, "/tumour/cna_calls_hgncsymbol.txt.gz")
  cn.matrix = read.delim(tumourFile, stringsAsFactors = F, row.names = 1)
  # Output number of samples with a cn alteration and number of samples used
  return(list(cna = (rowSums(cn.matrix != 0)) / ncol(cn.matrix), n = ncol(cn.matrix)))
})

# Number of samples with a cn alteration
cnProp.df = as.data.frame(t(ldply(lapply(cnProp.lst, function(x) return(x[[1]])))))
colnames(cnProp.df) = unlist(lapply(tumours, basename))

# Number of samples used
tumN.df = as.data.frame((ldply(lapply(cnProp.lst, function(x) return(x[[2]])))))
row.names(tumN.df) = unlist(lapply(tumours, basename))
write.table(x = tumN.df, file = "results/human/numberOfSamples_CNA_analysis.tab", quote = FALSE, sep = "\t")

# CNA heatmap
stress.properties.df <- read.delim("results/human/humanStressConsensus_Properties_Full_Table.tsv", stringsAsFactors = F)
# Remove those genes w/o entry in cna matrix
stress.properties.df.f = stress.properties.df[stress.properties.df$external_gene_name %in% row.names(cnProp.df),]

# HeLa
# Induced genes, order alphabetically
stressSel.ind = subset(stress.properties.df.f, salt_HeLa == "upregulated", select = c("external_gene_name"))
stressSel.ind = stressSel.ind[order(stressSel.ind$external_gene_name),]
# Repressed genes, order alphabetically
stressSel.rep = subset(stress.properties.df.f, salt_HeLa == "downregulated", select = c("external_gene_name"))
stressSel.rep = stressSel.rep[order(stressSel.rep$external_gene_name),]
# Unresponsive genes sample, order alphabetically
stressSel.unr = subset(stress.properties.df.f, salt_HeLa == "unresponsive", select = c("external_gene_name"))
set.seed(12345)
stressSel.unr = sample(stressSel.unr$external_gene_name, size = 450)
stressSel.unr = stressSel.unr[order(stressSel.unr)]

stressNew = c(stressSel.ind, stressSel.rep, stressSel.unr)

cnPropSel = cnProp.df[stressNew,]

# Converto to matrix
cnPropSel = as.matrix(cnPropSel)

rowSplit = c(rep("Upregulated", length(stressSel.ind)), 
             rep("Downregulated", length(stressSel.rep)),
             rep("Unresponsive", length(stressSel.unr)))

save(rowSplit, cnPropSel, file = "data/rdata/cna_tcga_HeLa.Rda")
write.table(x = cnPropSel, 
            file = "results/human/CNAproportion_TCGA_salt_HeLa.tab", 
            sep = "\t", 
            row.names = TRUE,
            quote = FALSE)
