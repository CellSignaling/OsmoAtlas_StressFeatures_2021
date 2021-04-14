# Function to compute Codon Pair Bias for yeast, as originally computed in:
# http://science.sciencemag.org/content/sci/suppl/2008/06/26/320.5884.1784.DC1/Coleman.SOM.pdf

# In this page you can download data in yeast for dicodon usage frequency that proves the computations are OK: 
# http://pare.sunycnse.com/cut/index.jsp
CodonPairBias <- function(genetable){
  library(seqinr)
  # Obtain for each gene: Table structure with dicodon, codon, aa and aa1aa2 count
  counts <- lapply(subset(genetable, !is.na(CDS))$CDS, function(cds){
    # Dicodon wordsize = 6 and an step of 3
    diCodons <- seqinr::count(s2c(tolower(cds)), wordsize = 6, by = 3)
    codons <- seqinr::count(s2c(tolower(cds)), wordsize = 3, by = 3)
    trans <- seqinr::translate(s2c(cds))
    aa <- seqinr::count(trans, wordsize = 1, by = 1,alphabet = a())
    diAa <- seqinr::count(trans, wordsize = 2, by = 1,alphabet = a())
    return(list(diCodons = diCodons, codons = codons, aa = aa, diAa = diAa))
  })
  # Data frame with rows = genes and columns = dicodons
  diCodonsDf  <- plyr::ldply(lapply(counts, with, diCodons))
  # Summarise in a single data frame with total counts for each dicodon
  diCodonsCount <- data.frame(diCodonCount = colSums(diCodonsDf))
  # Divide by the total to obtain frequency
  diCodonsCount$freq <- diCodonsCount$diCodonCount/sum(diCodonsCount$diCodonCount)
  
  # Codons
  codonsDf  <- plyr::ldply(lapply(counts, with, codons))
  codonsCount <- data.frame(codonCount = colSums(codonsDf))
  codonsCount$freq <- codonsCount$codonCount/sum(codonsCount$codonCount)
  
  # AA
  aaDf <- plyr::ldply(lapply(counts, with, aa))
  aaCount <- data.frame(aaCount = colSums(aaDf))
  aaCount$freq <- aaCount$aaCount/sum(aaCount$aaCount)
  
  # DiAA
  diAaDf <- plyr::ldply(lapply(counts, with, diAa))
  diAaCount <- data.frame(diAaCount = colSums(diAaDf))
  diAaCount$freq <- diAaCount$diAaCount/sum(diAaCount$diAaCount)
  
  codonsWoStop <- words()[!words() %in% c("tag", "taa", "tga")]
  
  # Obtain codon pairs combinations wo including stop codons
  codonsPairs <- expand.grid(codonsWoStop, codonsWoStop)
  codonsPairs$diCodons <- paste0(codonsPairs$Var1, codonsPairs$Var2)
  
  # For each codon pair obtain Codon Pair Score log(fxy/((fxfy*faa1aa2)/faa1*faa2)) 
  cpsLs <- lapply(codonsPairs$diCodons, function(pair){
    cd1 <- paste0(s2c(pair)[1:3], collapse = "")
    cd2 <- paste0(s2c(pair)[4:6], collapse = "")
    aa1 <- seqinr::translate(s2c(cd1))
    aa2 <- seqinr::translate(s2c(cd2))
    diAa <- paste0(seqinr::translate(s2c(pair)), collapse = "")
    denom <- (codonsCount[cd1,"freq"]*codonsCount[cd2,"freq"]*diAaCount[diAa, "freq"])/(aaCount[aa1, "freq"]*aaCount[aa2, "freq"])
    cps <- log(diCodonsCount[pair, "freq"]/denom)
  })
  
  # Generate a data frame with score and dicodon pair
  names(cpsLs) <- codonsPairs$diCodons
  cps <- plyr::ldply(cpsLs, data.frame)
  colnames(cps) <- c("dicodon", "cps")
  row.names(cps) <- cps$dicodon
  
  # Save data for reusing
  dir.create("data/rdata", showWarnings = F)
  save(cps, file = "data/rdata/cps.rda")
  
  cpb <- lapply(genetable$CDS, function(cds){
    if(is.na(cds)){
      print("No CDS")
      score = NA
      score
    } else {
      # Dicodon count for each CDS
      diCodons <- seqinr::count(s2c(tolower(cds)), wordsize = 6, by = 3)
      # Retain only dicodons that occur in the CDS in order to speed up the analysis
      diCodons <- diCodons[diCodons>0]
      # Retain only dicodons that appear in our cps table (that does not includes STOP codons)
      diCodons <- diCodons[names(diCodons) %in% row.names(cps)]

      score = sum(cps[names(diCodons),"cps"]*diCodons)/sum(diCodons)
      
      score
    }
  })
  
  names(cpb) <- row.names(genetable)
  save(cpb, file = "data/rdata/cpb.rda")
  
  cpbDf <- plyr::ldply(cpb, data.frame)
  colnames(cpbDf) <- c("gene", "cpb")
  row.names(cpbDf) <- cpbDf$gene
  
  return(cpbDf)
}
