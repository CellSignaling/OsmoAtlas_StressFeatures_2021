# Function to run the DESeq2 pipeline with desired settings (extended version with batch removal and consideration of batch groups)
# (lfc threshold, p-value threshold, betaPrior)
# that outputs PCA plot, volcano plot and results tsv file
## inputData: data frame with raw read counts
## samples: design matrix
## design: formula of the design (i.e ~condition or ~condition + genotype, etc)
## colGroups: name or names of columns in the design matrix that correspond to analyzed variables
## betaPrior: Default T for bulk RNA-seq experimets
## pAdjThreshold: threshold to consider a gene as signficant that it is used to optimal independent filtering
## lfcThreshold: threshold of lfc for statistical testing
## biomart: biomaRt object name, created before
## originalId: Format of id of genes in the raw counts
## destinationId: Format of if of destination after conversion with biomaRt
## genesInterest: String with names to be displayed in the graph, usually positive controls
## outputPrefix: Initial name of files produced by the analysis

DESeqQcResBatchExt <- function(inputData, 
                       samples,
                       design,
                       colGroups, 
                       variable,
                       betaPrior = T,
                       pAdjThreshold,
                       lfcThreshold,
                       biomart,
                       originalId,
                       destinationId,
                       genesInterest, 
                       outputPrefix
                       ){
  
  # Loading required libraries
  library(DESeq2)

  
  # Obtain second part of p-value threshold as a string (i.e 05 from 0.05)
  pThresholdStr <- unlist(strsplit(as.character(pAdjThreshold), ".", fixed = T))[2]
  
  # Extendend prefix for output file names
  outputPrefixExt <- paste0(outputPrefix, "_lfc", lfcThreshold, "_p", pThresholdStr)
  
  # Creating DESeq2 object 
  dds <- DESeqDataSetFromMatrix(countData = inputData, 
                                colData = samples, 
                                design = design)
  
  # Executing DESeq2 pipeline, with the option of set betaPrior to True (recomended for bulk RNA-seq, does not change much the result)
  dds <- DESeq(dds, betaPrior = betaPrior)

  # Regularized log transformation for different analysis (clustering, heatmaps, etc)
  rld <- rlogTransformation(dds)
  
  # Remove batch and plot
  rldBr = rld
  assay(rldBr) = limma::removeBatchEffect(assay(rldBr), rldBr$batch)
  
  # PCA
  # Extract rld counts and put them into a data frame
  rlogCounts <- as.data.frame(assay(rldBr))
  # Basal samples
  basal <- row.names(samples[samples[,colGroups] == "basal",])
  # Treatment samples
  treated <- row.names(samples[samples[,colGroups] == variable,])
  # Retain only considered samples
  rlogCountsS <- rlogCounts[,colnames(rlogCounts) %in% c(basal, treated)]
  # Only top500 by variance as in DESEq2 
  topVarIds = order(rowVars(as.matrix(rlogCountsS)), decreasing = T)[1:500]
  # Obtain PCs
  pcag <- prcomp(t(rlogCountsS[topVarIds,]))
  # Store PCs in a data frame
  scores <- as.data.frame(pcag$x)
  # Summary prpcomp
  info <- summary(pcag)
  infoVar <- as.data.frame(info$importance)
  # Varianze explained by each PC
  PC1Var <- round(infoVar$PC1[2] * 100)
  PC2Var <- round(infoVar$PC2[2] * 100)
  # Add group info
  sampData = as.data.frame(colData(rldBr))
  scores[,"group"] <- apply(sampData[,colGroups], 1, paste, collapse = ":")
  # Write data frame for plotting
  write.table(x = scores, file = paste0(outputPrefixExt, "_batchRemoved_pcagg_data.tsv"), sep = "\t", quote = F, row.names = T, col.names = NA)

  
  # Test for DE genes 
  res <- results(dds, 
                lfcThreshold = lfcThreshold, # Lfc testing changes the results (< 600 genes as compared to lfcThreshold = 0 and then selecting by lfc)
                altHypothesis = "greaterAbs", 
                alpha = pAdjThreshold)

  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  
  # Merge with normalizated count data
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort= FALSE)
  names(resdata)[1] <- "Gene"
  
  # Obtain gene symbol from ensembl id
  geneSymbols <- getBM(attributes = c(originalId, destinationId), 
        filters = originalId, 
        values = resdata$Gene, 
        mart = biomart)
  
  # Merge with results data frame
  resdataExt <- merge.data.frame(resdata, 
                   geneSymbols, 
                   by.x = "Gene", 
                   by.y = originalId)
  
  # Change name of ensembl id column from gene to Ensembl gene id
  colnames(resdataExt)[1] <- originalId
  
  # Reorder data frame for better understanding
  resdataExt <- cbind(resdataExt[, c(originalId, destinationId)], resdataExt[,!colnames(resdataExt) %in% c(originalId, destinationId)])
  
  # Write results
  write.table(resdataExt, 
              file= paste0(outputPrefixExt,"_results.tsv"), 
              quote = F, 
              sep = "\t", 
              row.names = F)
}
