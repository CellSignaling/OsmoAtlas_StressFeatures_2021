#######
#
#
# TCGA RNA-Seq processing and DESeq2 analysis
#
#######

# Loading required libraries/functions
library(dplyr)
library(purrr)
source("functions/TCGA_Rnaseq_analysis.R")

# TCGA folders
tcgaFiles = dir("data/original_data/tcga.rna/", full.names = T)

# List for results
tcgaRes = list()

# List for discarded tumors
discTum = list()

# List for number of samples analyzed
numSamp <- list()

# For each file do RNA-seq analysis
# Since it takes a lot of time, do it only if output does not exist

if(!file.exists("data/rdata/tcgaRes_TN_matched.rda")){

  for (cancer in tcgaFiles){
    # Extract type of cancer
    cancerType = unlist(lapply(strsplit(basename(cancer), split = ".", fixed = T), function(x) return(x[1])))
    # Print info message
    print(paste0("RNA-seq analysis of: ", cancerType))
    print("--- Reading and preprocessing")
    # Obtain cancertype RNA-Seq file
    cancerFile <- grep("patient_annotation", dir(cancer, full.names = TRUE), value = TRUE, invert = TRUE)
    # Read data
    rawData = read.delim(cancerFile, stringsAsFactors = F)
    # Remove scaled_estimate and transcript_id columns
    rawDataF = rawData[, !rawData[1,] %in% c("scaled_estimate", "transcript_id")]
    # Remove first row
    rawDataF = rawDataF[-1,]
    # Gene names as row names
    row.names(rawDataF) = rawDataF$Hybridization.REF
    # Remove gene names column
    rawDataF = rawDataF[,-1]
    # Transform data to numeric
    readyData = as.data.frame(sapply(rawDataF, as.numeric), row.names =  row.names(rawDataF))
    
    print("--- Tumor-normal matched samples")
    # Obtain tumor type information
    sampleType = as.integer(substr(colnames(readyData), 14, 15))
    # Select normal samples (11: Solid tissue normal)
    sampleTypeNormal = sampleType == 11
    
    # Sanity checking 1: No normal samples available
    if(all(!sampleTypeNormal)){
      
      warning(paste0("### ", cancerType," tumor type has no normal samples ###"))
      
      # Fill discarded tumors list
      discTum[cancerType] = "No normal samples"
      
    } else {
      
      # Select tumor matched samples (01: Primary solid tumor and same patient)
      patientIdNormal = substr(colnames(readyData)[sampleTypeNormal], 9,12)
      sampleTypeTumor = grepl(x = colnames(readyData), pattern = paste0(patientIdNormal, ".01", collapse = "|")) ## Note: This is buggy if patientIdNormal is empty
      
      # Sanity checking 2: No tumor samples available
      if(all(!sampleTypeTumor)){
        warning(paste0("### ", cancerType," tumor type has no normal-tumor matched samples ###"))
        
        # Fill discared tumors list
        discTum[cancerType] = "No normal-tumor matched samples"
        
      } else {
        # Filter unmatched normal
        patientIdTumor = substr(colnames(readyData)[sampleTypeTumor],9,12)
        sampleNormalMatched = grepl(x = colnames(readyData), pattern = paste0(patientIdTumor, ".11", collapse = "|"))
        
        # Filtered data frame
        readyData_normal_tumor = readyData[, sampleNormalMatched | sampleTypeTumor]
  
        nTotSamples = dim(readyData_normal_tumor)[2]
        print(paste0("--- Total number of collected samples: ", nTotSamples))
        # Store number of samples analyzed
        numSamp[cancerType] <- nTotSamples
        
        # Sanity checking  3: No correct matching normal-tumor
        if(all(table(substr(colnames(readyData_normal_tumor), 9, 12)) == 2)){
        # Do the analysis
          # Create sample annotation
          samples = data.frame(sample = colnames(readyData_normal_tumor), stringsAsFactors = F)
          samples <- samples %>%
            mutate(condition = ifelse(substr(sample, 14,15) == "11", yes = "normal", no = "tumor")) %>% # should work now but if is extended to other codes will fail
            mutate(patient = substr(sample, 9,12)) %>%
            tibble::column_to_rownames(var = "sample")
          samples$patient = factor(samples$patient)
          # Prepare dataset for Deseq2
          readyData_normal_tumor_counts = as.data.frame(sapply(readyData_normal_tumor, round), row.names = row.names(readyData_normal_tumor))
          # Check samples are in the same order as in data columns
          # all(row.names(samples)==colnames(readyData_normal_tumor_counts))
          # Execute DESEQ2 function
          outputPrefix = paste0("results/human/TCGA/", cancerType)
          print("--- Deseq2 analysis")
          resdata = minimalDeseq2(countData = readyData_normal_tumor_counts, 
                        sampAnno = samples, 
                        design = "~patient + condition",
                        parallel = F, 
                        pAdjThreshold = 0.05)
          # Plotting
          print("--- Plotting and formating output")
          # Obtain gene symbol from the row names
          resdata$geneSymbols <- unlist(lapply(strsplit(resdata$Gene, "|", fixed = T), function(x) return(x[1])))
          # NCBI ids
          resdata$ncbiIds <- unlist(lapply(strsplit(resdata$Gene, "|", fixed = T), function(x) return(x[2])))
          # Gene symbol if available if not ncbi
          resdata$identifier = ifelse(resdata$geneSymbols == "?", yes = resdata$ncbiIds, no = resdata$geneSymbols)
          # Look at "SLC35E2" annotation, 9906 is A anbd 728661 is B
          resdata[resdata$ncbiIds == 9906, c("Gene", "identifier")] = c("SLC35E2-A|9906", "SLC35E2-A")
          resdata[resdata$ncbiIds == 728661, c("Gene", "identifier")] = c("SLC35E2-B|728661", "SLC35E2-B")
          
          # Reorder data frame for better understanding
          resdataExt = resdata %>%
            select(Gene, geneSymbols, ncbiIds, identifier, everything())
          
  
          # Obtain data for next steps
          dataExt <- data.frame(gene = resdataExt[, "Gene"],
                                geneExt = resdataExt[, "identifier"],
                                padj = resdataExt$padj,
                                lpadj = -log10(resdataExt$padj),
                                lfc = resdataExt$log2FoldChange, 
                                stringsAsFactors = F)
          
          # Storing desired output in a list
          tcgaRes[[cancerType]] = dataExt[,c("geneExt", "padj", "lfc")]
        }
        else {
          warning("Tumor-normal matching error")
          
          # Fill discarded tumors list
          discTum[cancerType] <- "Tumor-normal matching error"
        }
      }
    }
  }
  
  
  
  ## Important, save results as soon as it finish
  save(tcgaRes, discTum, numSamp, file = "data/rdata/tcgaRes_TN_matched.rda")
}


# Load data
load("data/rdata/tcgaRes_TN_matched.rda", verbose = T)

# Add cancer name to lfc and pval colnames for later compatibility
for (i in names(tcgaRes)){
  colnames(tcgaRes[[i]]) = c("geneExt", paste0("padj_", i), paste0("lfc_", i))
  
}
# Join data in a data frame with shared keys (genes)
tcgaDf = tcgaRes %>%
  purrr::reduce(left_join, by = "geneExt")

# Intersect between experimental  data
## Read human stress genes tables
hsExp <- read.delim("results/human/humanStressConsensusTable_experimental_HeLa.tab", stringsAsFactors = F) 


# Prepair data for plotting
## Fold changes
foldChange = tcgaDf %>% 
  dplyr::select(geneExt, starts_with("lfc")) %>%
  tibble::column_to_rownames(var = "geneExt")

colnames(foldChange) = unlist(lapply(strsplit(colnames(foldChange), "_"), function(x) return(x[2])))

# Experimental genes: HeLa
# Filter those genes w/o correspondance in fold change table
stress.properties.df.f = hsExp[hsExp$associated_gene_name %in% row.names(foldChange),]
# Induced alphabetically
stressSel.ind = subset(stress.properties.df.f, salt_HeLa == "upregulated", select = c("associated_gene_name"))
stressSel.ind = stressSel.ind[order(stressSel.ind),]
# Repressed alphabetically
stressSel.rep = subset(stress.properties.df.f, salt_HeLa == "downregulated", select = c("associated_gene_name"))
stressSel.rep = stressSel.rep[order(stressSel.rep),]
# Unresponsive sample alphabetically
stressSel.unr = subset(stress.properties.df.f, salt_HeLa== "unresponsive", select = c("associated_gene_name"))
set.seed(12345)
stressSel.unr = sample(stressSel.unr$associated_gene_name, size = 450)
stressSel.unr = stressSel.unr[order(stressSel.unr)]


stressNew = c(stressSel.ind, stressSel.rep, stressSel.unr)

foldChangeSel = foldChange[stressNew,]


rowSplit = c(rep("Upregulated", length(stressSel.ind)), 
             rep("Downregulated", length(stressSel.rep)),
             rep("Unresponsive", length(stressSel.unr)))

save(rowSplit, foldChangeSel, file = "data/rdata/rnaseq_tcga_HeLa.Rda")
write.table(x = foldChangeSel, 
            file = "results/human/foldChange_TCGA_RNASeq_salt_HeLa.tab", 
            sep = "\t", 
            row.names = TRUE,
            quote = FALSE)
