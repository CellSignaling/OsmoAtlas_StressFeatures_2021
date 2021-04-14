# Script for computing gene expression variability from yeast single cell data (yscRNA-Seq)
###########################################################################################
# Loading required libraries
library(DESeq2)
library(scran)

# Data from Nadal-Ribelles, Islam, Wei, Latorre et al., 2019 (Supplementary table 3)
# Read raw data
rawSc = read.delim("data/original_data/rawExpression_YPD_NatMicro_w_erccs.tab", stringsAsFactors = F, sep = "\t")

# Remove HPH and KanMX (markers)
rawSc = rawSc[!rawSc$geneName %in% c("HPH", "KanMX"),]

# Deal with duplicated names
dIds = duplicated(rawSc$geneName)
rawSc$geneName[dIds] = paste(rawSc$geneName[dIds], "_2", sep = "")
row.names(rawSc) = rawSc$geneName

erids = grepl('ERCC', rawSc$geneName) # This works here, but careful with genes in human that start with ERCC


normDat.rcv2.dm.resVar <- function(expsfk, erids){
  # Obtain normalized read
  countsEndo <- expsfk[!erids,]
  sfEndo <- estimateSizeFactorsForMatrix(expsfk[!erids,])
  nCountsEndo <- t( t(countsEndo) / sfEndo )
  row.names(nCountsEndo) <- rawSc$geneName[!erids]
  
  # Obtain normalized read for ERCCs
  countsERCC <- expsfk[erids,]
  sfERCC <- estimateSizeFactorsForMatrix( expsfk[erids,] )
  nCountsERCC <- t( t(countsERCC) / sfERCC )
  row.names(nCountsERCC) <- rawSc$geneName[erids]
  
  # Join counts
  nCounts <- rbind(nCountsEndo, nCountsERCC)
  
  # Compute mean and variance for each gene
  means <- rowMeans(nCounts)
  vars <- rowVars(nCounts)
  
  # Create a data frame with all the information and the computed cv2
  ysc.gene.summary = data.frame(name = row.names(nCounts), Mean = means, Var = vars)
  ysc.gene.summary$cv2 <- ysc.gene.summary$Var/(ysc.gene.summary$Mean)^2
  
  # Obtain DM
  ysc.gene.summary$DM <- DM(mean = ysc.gene.summary$Mean, cv2 = ysc.gene.summary$cv2)
  
  # scTransform computation (input is UMI) #
  # Obtain output from vst
  vst_out <- sctransform::vst(as.matrix(expsfk), latent_var = c('log_umi'), return_cell_attr = T, verbosity = T)
  ysc.gene.summary$vst_residual_variance <- vst_out$gene_attr[match(ysc.gene.summary$name, row.names(vst_out$gene_attr)), "residual_variance"]
  
  # Table for output
  ysc.gene.out <- ysc.gene.summary[,c("name", "DM", "vst_residual_variance")]
  
  # Write output
  write.table(x = ysc.gene.out, 
              file = "results/S_cerevisiae/expsVariability_DM_resVar.tab", 
              sep = "\t", 
              quote = F, 
              row.names = F)
}


normDat.rcv2.dm.resVar(expsfk = rawSc[, -c(1,2)], erids = erids)
