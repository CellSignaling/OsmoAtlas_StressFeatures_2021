###################
#
# Script for gathering features/properties of human genes
# 
###################

##### Gather properties of human genes #####

# Make it as much similar as possible to the yeast project

# Load required libraries/functions
library(biomaRt)
library(stringr)
source("functions/CodonPairBias.R")
library(qusage) # read.gmt
library(dplyr)
library(seqinr)

# Load ENSEMBL
hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://may2017.archive.ensembl.org")

# Read identifiers for compatibility
humanIds = read.delim("data/derived_data/hs_GRCh38.p13_geneNames_singleTranscript.tab", stringsAsFactors = F)

# Read data
geneswUtrs = read.delim("data/derived_data/hs_GRCh38.p13_w_UTRs_sequences_r_format.tab", stringsAsFactors = F, row.names = 1)

# Property 1: Length
# lapply gives the result to a list
# with do.call, as.data.frame and cbind
# ordering as a data frame, specifying
# row.names from geneswUTRs
lengthGeneswUtrs <- as.data.frame(do.call(cbind, lapply(geneswUtrs, nchar)), 
                                  row.names = row.names(geneswUtrs))
colnames(lengthGeneswUtrs) <- paste0("length_",colnames(lengthGeneswUtrs))

# Property 2: GC content comparison
gcGeneswUtrs <- as.data.frame(do.call(cbind, lapply(geneswUtrs, function(x) (str_count(x, "G") + str_count(x, "C"))/ nchar(x))), 
                              row.names = row.names(geneswUtrs))
colnames(gcGeneswUtrs) <- paste0("gc_",colnames(gcGeneswUtrs))


## Property 3: Codon Pair Bias
# http://science.sciencemag.org/content/sci/suppl/2008/06/26/320.5884.1784.DC1/Coleman.SOM.pdf
cpbGeneswUtrs <- CodonPairBias(genetable = geneswUtrs)
cpbGeneswUtrs <- data.frame(cpbGeneswUtrs, row.names = 1, stringsAsFactors = F)
save(cpbGeneswUtrs, file = "data/rdata/cpb_human.rda")
load("data/rdata/cpb_human.rda", verbose = TRUE)


## Property 4: Number of TF
# ENCODE based on Chip-seq (n = 22819)
encodeTf = read.gmt("data/original_data/harmonizome_ENCODE_TF_attribute_set_library_crisp.gmt")
nTfGeneswUtrs = plyr::ldply(lapply(encodeTf, length), as.data.frame)
colnames(nTfGeneswUtrs) = c("hsapiens_associated_gene_name", "nTf")

# Transform to ensembl ids
nTfGeneswUtrs$ensembl_gene_id = humanIds[match(nTfGeneswUtrs$hsapiens_associated_gene_name, humanIds$external_gene_name), "ensembl_gene_id"]
# Remove those without an Ensembl id entry
nTfGeneswUtrs = nTfGeneswUtrs[!is.na(nTfGeneswUtrs$ensembl_gene_id),]
# Put ensembl as row names
row.names(nTfGeneswUtrs) = nTfGeneswUtrs$ensembl_gene_id
# Remove external name 
nTfGeneswUtrs = nTfGeneswUtrs %>%
  dplyr::select(-c(hsapiens_associated_gene_name, ensembl_gene_id))


## Property 5: Halflife
## Tani dataset (BRIC-seq in HeLa)
# Modifed dataset: >24 to 24 and N.D. to NA
halfLife.tani = read.delim("data/original_data/halfLife_hr_Tani_et_al_2012_Gen_research.txt", stringsAsFactors = F)
# Retain only one identifier
halfLife.tani$RepName = unlist(lapply(strsplit(halfLife.tani$RepName, ","), function(x) return(x[1])))
RefseqToEns = getBM(attributes = c("refseq_mrna","ensembl_gene_id"), filters = "refseq_mrna", values =  halfLife.tani$RepName,mart = hs)
# Eliminate duplicated Refseq
RefseqToEns = RefseqToEns[!duplicated(RefseqToEns$refseq_mrna),]
# Eliminate duplicated ensembl
RefseqToEns = RefseqToEns[!duplicated(RefseqToEns$ensembl_gene_id), ]
# Map back to the halflife data frame
halfLife.tani$ensembl_gene_id = RefseqToEns[match(halfLife.tani$RepName, RefseqToEns$refseq_mrna), "ensembl_gene_id"]
# Filter NAs
halfLife.tani = halfLife.tani[!is.na(halfLife.tani$ensembl_gene_id),]
halfLife.tani = halfLife.tani[!is.na(halfLife.tani$halfLife_tani),]

halfLifeGeneswUtrs = halfLife.tani
row.names(halfLifeGeneswUtrs) = halfLifeGeneswUtrs$ensembl_gene_id

halfLifeGeneswUtrs = halfLifeGeneswUtrs %>%
  dplyr::select(-c(ensembl_gene_id, RepName))

colnames(halfLifeGeneswUtrs) <- "halfLife"
# Property 6: Distance to Median
## HeLa cells
#### In this case there are no ERCCs, if you do  grep "ERRC" you find some genes but careful
#### because they are Excision Repair proteins
scHela_1 = read.delim("data/original_data/GSM3713084_HeLa_1.txt", stringsAsFactors = F)
scHela_2 = read.delim("data/original_data/GSM3713085_HeLa_2.txt", stringsAsFactors = F)
scHela_3 = read.delim("data/original_data/GSM3713086_HeLa_3.txt", stringsAsFactors = F)

# Duplicated row names in scHela_3 are 1-Mar and 2-Mar, remove them
scHela_3 = scHela_3[!scHela_3$X %in% c("1-Mar", "2-Mar"),]

commonIds = Reduce(intersect, list(scHela_1$X, scHela_2$X, scHela_3$X))

tmpH1 = scHela_1[scHela_1$X %in% commonIds, ]
tmpH2 = scHela_2[scHela_2$X %in% commonIds, ]
tmpH3 = scHela_3[scHela_3$X %in% commonIds, ]

scHela = Reduce(left_join, list(tmpH1, tmpH2, tmpH3))
scHela = tibble::column_to_rownames(scHela, var = "X")

# Filter cells using criteria from the original paper (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0225466)
## i.e filter cells with less than 2000 detected genes
nDetectGenes = apply(scHela, 2,function(cell){sum(cell != 0)})
scHela.f = scHela[, !nDetectGenes <= 2000] # Number of cells after filtering is the same as stated on the paper

## Filter genes with no detected expression in the cells left
gfIds = rowSums(scHela.f) == 0
scHela.f = scHela.f[!gfIds,]


# Not ERCCs available
# Obtain normalized data
countsEndo = scHela.f
sfEndo = DESeq2::estimateSizeFactorsForMatrix(countsEndo)
nCounts = t( t(countsEndo) / sfEndo )

# Compute mean and variance for each gene
means = rowMeans(nCounts)
vars = matrixStats::rowVars(nCounts)

# Create a data frame with all the information and the computed cv2
gene.summary = as.data.frame(cbind(means, vars))
colnames(gene.summary) = c("Mean", "Var")
gene.summary$cv2 = gene.summary$Var/((gene.summary$Mean)^2)

# supplementary Note 7 - On "maximal" noise from Brenecke et al., 2013
## Basically: For  a  given  gene  with average normalized read count μ, the maximum CV! is reached when 
# all but one of the mcells  have  zero  counts  and  one  sample  has  a  normalized  count  of mμ,  
# resulting  in  a variance  of mμ!  and  hence  a  CV!  of m

# Briefly, a median-based trend is fitted to the log-transformed cv2 against the log-transformed mean using runmed
gene.summary$DM <- scran::DM(gene.summary$Mean, gene.summary$cv2)

gene.summary$ensembl_gene_id = humanIds[match(row.names(gene.summary), humanIds$external_gene_name),"ensembl_gene_id"]

dmGeneswUtrs = gene.summary[, c("ensembl_gene_id","DM")]
# Remove those without an Ensembl id entry
dmGeneswUtrs = dmGeneswUtrs[!is.na(dmGeneswUtrs$ensembl_gene_id),]
row.names(dmGeneswUtrs) = dmGeneswUtrs$ensembl_gene_id
# Remove ensembl gene id from column 
dmGeneswUtrs = dmGeneswUtrs %>%
  dplyr::select(-c(ensembl_gene_id))


# Property 7: Broad conservation: Human data
# Load Human InParanoid library
library("hom.Hs.inp.db")
# Making a connection to the database
# Why not to use preloaded InParanoid objects? Because they only consider ortholog-seed pairs with a 100% bootstrap score
# In Koch et al., 2012 they consider apparently all ortholog-seed pairs so we do so
# Due to that, we need to connect to the database (version 7.0, 100 species, now is default)
# Load library for connecting with databases
library(DBI)
# Create connection with the database
mycon <- hom.Hs.inp_dbconn()
# List all the table available in the database (i.e the species plus some metadata)
species = dbListTables(mycon)
# Filter metadata
species = species[!species %in% c("map_counts", "map_metadata", "metadata")]
# In here we are not filtering yeast species since there is not yeast conservation approach
# Generate a query for all the species
# The query is selecting all the elements in a a given (species) table
# This consist about human and X species genes
# If we select those genes that are human and that have bootstrap score, this mean those genes are all seeds in an ortholog-seed pair
sqlList = paste0("SELECT * FROM ", species)

genesWithOrthologBySp = lapply(sqlList, function(sql){
  # Get from the database the species table
  dataOut = dbGetQuery(mycon, sql)
  # Retain only human seeds of the ortholog pairs 
  tmp = subset(dataOut, seed_status != "" & species == "HOMSA")
  # Return the Inparanoid ID (ENSEMBL protein ID)
  return(tmp$inp_id)
})

# Count total number of ocurrences of a given gene in the 99 species (100 counting on Hs).
# In essence, we are counting the number of species in which an human gene has an ortholog
broadConservationGeneswUtrs = plyr::ldply(table(unlist(genesWithOrthologBySp)))
colnames(broadConservationGeneswUtrs) = c("ensembl_peptide_id", "Broad_conservation")

# We don't know if the rest of genes that are not in any cluster
# they are not because they have no orthologs in any species or because they have not
# been checked by the InParanoid algorithm, then we will assign NA or zero depending.
# However, this could be checked by looking at the FASTA sequences used in the version 7 of the database
# Data downloaded from: http://inparanoid.sbc.su.se/download/old_versions/data_7.0/sequences.tgz
inpV7Hs <- read.fasta(file = "data/original_data/H.sapiens_inparanoid_v7_processed.fa")
inpV7HsIds <- names(inpV7Hs)
# We can assign zero to the proteins that were considered in Inparanoid v7 but don't have orthologs
# NA to the rest that does not appear here
zeroConsIds <- inpV7HsIds[!inpV7HsIds %in% broadConservationGeneswUtrs$ensembl_peptide_id]
zeroCons <- data.frame(ensembl_peptide_id = zeroConsIds, Broad_conservation = 0, stringsAsFactors = FALSE)

broadConservationGeneswUtrs <- rbind(broadConservationGeneswUtrs, zeroCons)

ensProt_Gene = getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"), filters = "ensembl_peptide_id", values = broadConservationGeneswUtrs$ensembl_peptide_id, mart = hs)
# Remove duplicates (i.e several proteins mapping to the same gene)
ensProt_Gene = ensProt_Gene[!duplicated(ensProt_Gene$ensembl_gene_id),]

# Prepair data with the desired format
broadConservationGeneswUtrs$ensembl_gene_id = ensProt_Gene[match(broadConservationGeneswUtrs$ensembl_peptide_id, ensProt_Gene$ensembl_peptide_id), "ensembl_gene_id"]
broadConservationGeneswUtrs = broadConservationGeneswUtrs[!is.na(broadConservationGeneswUtrs$ensembl_gene_id),]

row.names(broadConservationGeneswUtrs) = broadConservationGeneswUtrs$ensembl_gene_id
broadConservationGeneswUtrs = broadConservationGeneswUtrs %>%
  dplyr::select(-c(ensembl_gene_id, ensembl_peptide_id))


# Property 8: Disorder --> Percentage of disordered residues in a protein
# Using D2P2 database and settings from van der lee et al 2014., chemical reviews
# Residues in each protein are defined as disordered when there is a consensus between >75% of the predictors in the D2P2d at a base at that position. 
# The set of human genes was taken from Ensembl release63 (Disopred page), 
# and the representative protein coded for by the longest transcript was used in each case
# Download the IDS used by D2P2 from their webpage: http://d2p2.pro/downloads/genomes.protein.gz
# awk '$1 == "hs"' data/original_data/d2p2_genomes.protein  > data/original_data/hs_d2p2_enspIds_ensmblv63.txt
d2p2Ids = read.delim("data/original_data/hs_d2p2_enspIds_ensmblv63.txt", stringsAsFactors = F, col.names = c("org","ensembl_peptide_id", "n", "desc"))

# Obtain gene and transcript identifier for later use
d2p2Ids$ensembl_gene_id = unlist(lapply(strsplit(unlist(lapply(strsplit(d2p2Ids$desc, " "), function(x) return(tail(x, 2)[1]))), ":"), function(y) return(y[2])))
d2p2Ids$ensembl_transcript_id = unlist(lapply(strsplit(unlist(lapply(strsplit(d2p2Ids$desc, " "), function(x) return(tail(x, 1)))), ":"), function(y) return(y[2])))

# Since it takes a lot of time, do it only if output does not exist

if(!file.exists("data/rdata/disorder_d2p2.Rda")){
  # For each protein id:
  library('rjson')
  library(plyr)
  library(pbapply)
  seqIds = paste0('http://d2p2.pro/api/seqid/["', d2p2Ids$ensembl_peptide_id,'"]')
  seqDisorder = pblapply(seqIds, cl = 5, function(indSeq){  # Use pblapply parallel
    print(indSeq)
    # Obtain JSON format from d2p2 webpage
    response = fromJSON(readLines(indSeq, warn = FALSE))
    # Return NA if there is no info for a given protein id
    if(length(response[[1]]) == 0){
      return(c(NA,NA))
    } else {
      # Obtain the ranges in which the protein is disorderd
      consRanges = response[[1]][[1]][[3]]$disorder$consranges
      nDisRes = sum(unlist(lapply(consRanges, function(y){
        (as.numeric(y[2]) - as.numeric(y[1])) + 1
      })))
      # Obtain the total number of residues of a sequence
      totRes = length(response[[1]][[1]][[3]]$disorder$consensus)
      return(c(nDisRes, totRes))
    }
  })
  
  # Save as it takes approx 5h to compute
  save(seqDisorder, file = "data/rdata/disorder_d2p2.Rda")
}
load("data/rdata/disorder_d2p2.Rda", verbose = T)

d2p2Ids[,c("nDisRes", "totRes")] <- plyr::ldply(seqDisorder)

# Remove gene that has 0 length due to a bug in d2p2 database (i.e has no information about disorder in the downloadable file)
d2p2Ids = d2p2Ids[d2p2Ids$ensembl_peptide_id != "ENSP00000447029",]

# Compute proportion of disordered residues per protein
d2p2Ids$Disorder = d2p2Ids$nDisRes/d2p2Ids$totRes


# Compare with Figure 2 from https://pubs.acs.org/doi/pdf/10.1021/cr400525m; hist(d2p2Ids$Disorder*100, 11). Very similar

# Retain only transcripts that are in the reference data frame
disorderGeneswUtrs = d2p2Ids[d2p2Ids$ensembl_transcript_id %in% humanIds$ensembl_transcript_id,]
# Remove duplicated genes entries (randomly)
disorderGeneswUtrs = disorderGeneswUtrs[!duplicated(disorderGeneswUtrs$ensembl_gene_id),]
row.names(disorderGeneswUtrs) = disorderGeneswUtrs$ensembl_gene_id
disorderGeneswUtrs = disorderGeneswUtrs %>%
  dplyr::select(c(Disorder))


# Property 9: PPI degree
# Download BIOGRID from: https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-3.5.181/BIOGRID-ORGANISM-3.5.181.tab2.zip
hs_biogrid = read.delim("data/original_data/BIOGRID-ORGANISM-Homo_sapiens-3.5.181.tab2.txt", stringsAsFactors = F)
# Physical interactions
hs_biogrid = subset(hs_biogrid, Experimental.System.Type == "physical")
# Filter out "Affinity Capture-Luminescence" and "Reconstituted Complex" as Koch et al., 2012
## Careful: In Koch et al., 2012 they consider affinity-RNA. That is a protein-RNA interaction. Take it for consistency but take into account
hs_biogrid = hs_biogrid[!hs_biogrid$Experimental.System %in% c("Affinity Capture-Luminescence", "Reconstituted Complex"),]
# Filter out interactions with other organisms (i.e retain only human)
hs_biogrid = hs_biogrid[hs_biogrid$Organism.Interactor.A == 9606  & hs_biogrid$Organism.Interactor.B == 9606,]
# Filter out innecesary information
hs_biogrid = subset(hs_biogrid,,c( "X.BioGRID.Interaction.ID","Entrez.Gene.Interactor.A", "Entrez.Gene.Interactor.B","Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B"))
# Find interactors for each gene
geneInteractorsLst = with(hs_biogrid, split(Entrez.Gene.Interactor.B, Entrez.Gene.Interactor.A))
geneInteractorsLst2 = with(hs_biogrid, split(Entrez.Gene.Interactor.A, Entrez.Gene.Interactor.B))
# In Biogrid there is Bait and prey (directionality) in the interactions. The same protein could be bait and pray. So,
# we consider all the interactors of a protein as a bait and as a pray for the final number of interactors
keys <- unique(c(names(geneInteractorsLst), names(geneInteractorsLst2)))
geneInteractorsTot = setNames(mapply(c, geneInteractorsLst[keys], geneInteractorsLst2[keys]), keys)
# Count number of interactors per gene
geneInteractorsCount = lapply(geneInteractorsTot, function(gene) length(unique(gene)))
# Generate a data frame
ppiDegreeDf = plyr::ldply(geneInteractorsCount)
colnames(ppiDegreeDf) = c("entrez_gene_id", "PPI_degree")
# Convert from NCBI to ENSEMBL data frame ready
entrezToEnsembl = getBM(attributes = c("entrezgene", "ensembl_gene_id"), 
                        filters = c("entrezgene"), 
                        values = ppiDegreeDf$entrez_gene_id, 
                        mart = hs)

# Remove those genes that are not included in our primary list of ensembl identifiers with just one transcript isoform selected
entrezToEnsembl.f = entrezToEnsembl[entrezToEnsembl$ensembl_gene_id %in% humanIds$ensembl_gene_id,]
# Remove duplicates (randomly)
entrezToEnsembl.f = entrezToEnsembl.f[!duplicated(entrezToEnsembl.f$entrezgene),]
entrezToEnsembl.f = entrezToEnsembl.f[!duplicated(entrezToEnsembl.f$ensembl_gene_id),]
# Match ensembl with entrez in the ppi dataframe
ppiDegreeDf$ensembl_gene_id = entrezToEnsembl.f[match(ppiDegreeDf$entrez_gene_id, entrezToEnsembl.f$entrezgene), "ensembl_gene_id"]
# Prepare final data frame
ppiDegreeGeneswUtrs = ppiDegreeDf[!is.na(ppiDegreeDf$ensembl_gene_id ),]
row.names(ppiDegreeGeneswUtrs) = ppiDegreeGeneswUtrs$ensembl_gene_id
ppiDegreeGeneswUtrs = ppiDegreeGeneswUtrs %>%
  dplyr::select(-c("entrez_gene_id", "ensembl_gene_id"))

# Gathering properties in a single table
# Doing the gather of the data frames in an automatic way
properties.dataframes = grep("GeneswUtrs", ls(), value = T)
properties.list <- lapply(properties.dataframes, get)
properties.list <- lapply(properties.list, function(x) data.frame(x, ensembl_gene_id = row.names(x), stringsAsFactors = F))
names(properties.list) <- properties.dataframes
save(properties.list, file = "data/rdata/allProperties_hsapiens_unformated_list.Rda")
load("data/rdata/allProperties_hsapiens_unformated_list.Rda", verbose = T)
# Creating an empty data frame with the identifieres of the 56487 reference genes annotation (geneswUtrs)
# In this way we ensure to have a consistent annotation with the genes that match
properties.df <- data.frame("ensembl_gene_id" = humanIds$ensembl_gene_id, "external_gene_name" = humanIds$external_gene_name, stringsAsFactors = F)
print("Merging the following properties into a single data frame:")
for (i in names(properties.list)){
  property.name <- unlist(strsplit(i, "GeneswUtrs"))
  print(property.name)
  property <- properties.list[[i]]
  properties.df <- merge(properties.df, property, by = "ensembl_gene_id", all.x = T)
}
# Count number of complete observations for each property
nObservations <- plyr::ldply(colSums(!sapply(properties.df, is.na)))
colnames(nObservations) <- c("Property", "Observations")

# Read and join gene tables
humanGenesExp = read.delim("results/human/humanStressConsensusTable_experimental_HeLa.tab", stringsAsFactors = F)

humanGenesExp = humanGenesExp %>%
  dplyr::select(-c("associated_gene_name"))
colnames(humanGenesExp)[1] = "ensembl_gene_id"

stressConsensusProperties = merge(properties.df, humanGenesExp, by = "ensembl_gene_id", all.x = T)


# Order data frame
stressConsensusProperties = stressConsensusProperties %>%
  dplyr::select("ensembl_gene_id", "external_gene_name", "salt_HeLa", everything())

# Write output
write.table(x = stressConsensusProperties, file = "results/human/humanStressConsensus_Properties_Full_Table.tsv", sep = "\t", quote = F, row.names = F)
