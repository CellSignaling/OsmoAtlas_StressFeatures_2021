minimalDeseq2 <- function(countData, sampAnno, design, pAdjThreshold = 0.05, parallel = T){
  library(DESeq2)
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = sampAnno,
                                design = as.formula(design))
  # Executing DESeq2 pipeline, with the option of set betaPrior to True (recomended for bulk RNA-seq)
  dds <- DESeq(dds, betaPrior = T, parallel = parallel)

  # Obtain results with a given cutoff
  res = results(dds, alpha = pAdjThreshold, parallel = parallel)
  
  # Merge with normalizated count data
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort= FALSE)
  names(resdata)[1] <- "Gene"
  
  # Return table as output from the function
  return(resdata)
}
