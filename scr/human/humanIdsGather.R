###################
#
# Script for obtainzing human gene ids (1 gene - 1 transcript)
# 
###################


# Load required libraries

library(biomaRt)
library(dplyr)
hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://may2017.archive.ensembl.org")

# Download gene IDS from biomart
humanIds = getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "transcript_appris"), 
                 mart = hs)

# Select those that only have a transcript (39792 out of 63898)
txByGene = with(humanIds, split(ensembl_transcript_id, ensembl_gene_id))
singIsoTxs = names((which(lapply(txByGene, length) == 1)))
singHumanIds = humanIds[humanIds$ensembl_gene_id %in% singIsoTxs,]

# Filter then from humanIds to do the selection among the rest
humanIdsF = humanIds[!humanIds$ensembl_gene_id %in% singIsoTxs, ]

# Selection with APPRIS
# First: principal1
pral1 = subset(humanIdsF, transcript_appris == "principal1")
pral1 = pral1[!duplicated(pral1$ensembl_gene_id),]
# remove from humanIdsF
humanIdsFF = humanIdsF[!humanIdsF$ensembl_gene_id %in% pral1$ensembl_gene_id,]
# Second: principal2
pral2 = subset(humanIdsFF, transcript_appris == "principal2")
pral2 = pral2[!duplicated(pral2$ensembl_gene_id),]
# remove from humanIdsF
humanIdsFFF = humanIdsFF[!humanIdsFF$ensembl_gene_id %in% pral2$ensembl_gene_id,]
# Third: principal3
pral3 = subset(humanIdsFFF, transcript_appris == "principal3")
pral3 = pral3[!duplicated(pral3$ensembl_gene_id),]
# remove from humanIdsF
humanIdsFFFF = humanIdsFFF[!humanIdsFFF$ensembl_gene_id %in% pral3$ensembl_gene_id,]
# Fourth: principal4
pral4 = subset(humanIdsFFFF, transcript_appris == "principal4")
pral4 = pral4[!duplicated(pral4$ensembl_gene_id),]
# remove from humanIdsF
humanIdsFFFFF = humanIdsFFFF[!humanIdsFFFF$ensembl_gene_id %in% pral4$ensembl_gene_id,]

# Join data frames
humanGenesSingleTx = rbind(singHumanIds, pral1, pral2, pral3, pral4) 

# Get sequences of transcripts (unique transcripts)
utr3 <- getSequence(id=humanGenesSingleTx$ensembl_transcript_id, type="ensembl_transcript_id", seqType="3utr", mart=hs)
utr5 = getSequence(id=humanGenesSingleTx$ensembl_transcript_id, type="ensembl_transcript_id", seqType="5utr", mart=hs)
cds = getSequence(id=humanGenesSingleTx$ensembl_transcript_id, type="ensembl_transcript_id", seqType="coding", mart=hs)

save(utr3, utr5, cds, file = "data/rdata/sequences.Rda")

load("data/rdata/sequences.Rda", verbose = T)

# Join UTRs and CDS
humanSequences = Reduce(function(x,y) merge(x,y, by = "ensembl_transcript_id"), list(utr5, cds, utr3))
# Extra column with gene identifier and gene symbol
genesWUtrs = merge(humanGenesSingleTx[,1:3], humanSequences, by = "ensembl_transcript_id")

# For extra match with yeast work, substitute Sequence unavailable by NA and put only gene names
genesWUtrs[genesWUtrs == "Sequence unavailable"] = NA
write.table(genesWUtrs[,1:3], file = "data/derived_data/hs_GRCh38.p13_geneNames_singleTranscript.tab", sep = "\t", quote = F, row.names = F)

genesWUtrsFinal = genesWUtrs %>% dplyr::select(-c(ensembl_transcript_id, external_gene_name))
colnames(genesWUtrsFinal) = c("ensembl_gene_id", "UTR5", "CDS", "UTR3")
write.table(genesWUtrsFinal, file = "data/derived_data/hs_GRCh38.p13_w_UTRs_sequences_r_format.tab", sep = "\t", quote = F, row.names = F)
