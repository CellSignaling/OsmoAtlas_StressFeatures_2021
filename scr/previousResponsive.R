# Process stress-responsive genes from Gasch et al., 2000 and Nadal-Ribelles et al., 2014

## Gasch et al., 2000

outGasch <- "data/original_data/responsive_Gasch.txt"

# Download data
if (!file.exists(outGasch)){
  system2(command = "wget", args = c("http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt", 
                                     "-O", 
                                     outGasch))  
}

resGasch <- read.delim(outGasch, 
                       stringsAsFactors = F, 
                       check.names = F)

# 1 and 2 are UID and NAME, respectively
# The data contained in all files represents the normalized, 
# background-corrected log2 values of the Red/Green ratios measured on the DNA microarrays.
# select only sorbitol samples

osmoResGasch <- resGasch[,c(1, 2,grep("sorbitol", colnames(resGasch)))]

# select sorbitol 15min and write the output
osmoFc = subset(osmoResGasch, ,c(UID,`1M sorbitol - 15 min`))

write.table(osmoFc, 
            file = "data/derived_data/1M_Sorbitol_15min_FC_Gasch.txt", 
            quote = F, 
            row.names = F, 
            col.names = T, 
            sep = "\t")


## Tiling arrays, from Mariona directly (Nadal-Ribelles et al., 2014). The data here is in FC. 

resTilings <- read.delim("data/original_data/responsive_tilings.txt", 
                          stringsAsFactors = F, 
                          check.names = F)

# Only ORFs and from SteinmetzLab
resTilingsF <- subset(resTilings, source=="SteinmetzLab" & type == "ORF-T")

# Write full FC for 0.4M, 15m
tilingsFc = subset(resTilingsF, , c(name,`wt 15' 0.4M - wt 0`))
write.table(tilingsFc, 
            "data/derived_data/0.4M_NaCl_FC_15min_responsive_Tilings.txt", 
            quote = F, 
            col.names = T, 
            row.names = F, 
            sep = "\t")