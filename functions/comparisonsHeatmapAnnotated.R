# Function that performs pairwise wilcoxon test
# comparing stress responsive groups against unresponsive genes
# for all the properties. Then, summarizes this into a heatmap with an annotation
# that shows if a comparison is significative or not and 
# the direction of change. The order of the properties in the heatmap is alphabetic
# ARGUMENTS
# propertiesStressData = data frame with the properties and the stress genes classification
# properties = Properties to be compared 
# columnDep = Condition to be compared in terms of stress responsive genes
# random = wether to choose a given number of random genes in the unresponsive group
# sourceMolTable = Table with info about source and molecule of the properties
# outWidth = width of the output plot
# outHeight = height of the output plot
# outTable = name of the table with the results of the analysis
comparisonsHeatmapAnnotated <- function(propertiesStressData, 
                                        properties, 
                                        columnDep, 
                                        random = 0,
                                        outWidth, 
                                        outHeight,
                                        sourceMolTable,
                                        outTable){
# Loading required libraries
library(wesanderson)
library(pheatmap)  

if(random!=0){
  set.seed(12345)
  unresponsiveSample <- sample(which(propertiesStressData[, columnDep] == "unresponsive"), size = random)
  responsive <- which(propertiesStressData[,columnDep] != "unresponsive")
  propertiesStressData <- propertiesStressData[c(unresponsiveSample, responsive), ]
  outputPrefix <- paste0("_", columnDep,"_random_n_", random)
} else if (random == 0){
  outputPrefix <- paste0("_", columnDep)
}

## Obtain p-values of a Wilcoxon test
print("Computing p-values of Wilcoxon test for each property")

pvals.list <- lapply(properties, function(prop){
  property.group.df <- propertiesStressData[,c(columnDep,prop)]
  wt <- pairwise.wilcox.test(property.group.df[,prop], property.group.df[,columnDep], p.adjust.method = "none")
  wt.pvals <- c(wt$p.value["unresponsive", "downregulated"], wt$p.value["upregulated", "unresponsive"])
  names(wt.pvals) <- c("downregulated", "upregulated")
  return(wt.pvals)
})
names(pvals.list) <- properties
pvals.df <- plyr::ldply(pvals.list, .id = "Property")

## p-values correction
pvals.adj.df <- as.data.frame(matrix(p.adjust(as.vector(as.matrix(pvals.df[, c("downregulated", "upregulated")])), method='fdr'),ncol=2))
pvals.adj.df <- cbind(pvals.df$Property, pvals.adj.df)
colnames(pvals.adj.df) <- c("Property", "downregulated", "upregulated")


## Obtain comparisons of medians
print("Computing and comparing medians for each property")

medians.comp.list <- lapply(properties, function(prop){
  property.group.df <- propertiesStressData[,c(columnDep,prop)]
  property.group.list <- split(property.group.df[,prop], property.group.df[,columnDep])
  medians.list <- lapply(property.group.list, function(x) median(x, na.rm = T))
  medians.comp <- ifelse(unlist(medians.list) > medians.list$unresponsive, 
                         yes = "Higher", no = ifelse(unlist(medians.list) == medians.list$unresponsive, 
                                                     yes = "Equal", no = "Lower"))
  medians.comp <- medians.comp[grep("unresponsive", names(medians.comp), invert = T)]
  return(medians.comp)
})
names(medians.comp.list) <- properties
medians.comp.df <- plyr::ldply(medians.comp.list, .id = "Property")

# Compute median values for output table
medians.values.list <- lapply(properties, function(prop){
  property.group.df <- propertiesStressData[,c(columnDep,prop)]
  property.group.list <- split(property.group.df[,prop], property.group.df[,columnDep])
  medians.list <- lapply(property.group.list, function(x) round(median(x, na.rm = T), digits = 3))
  return(unlist(medians.list))
})
names(medians.values.list) <- properties
medians.values.df <- plyr::ldply(medians.values.list, .id = "Property")


medians.out.df <- merge(medians.values.df, pvals.adj.df, by = "Property")
colnames(medians.out.df) <- c("Property", "median_Downregulated", "median_Unresponsive", "median_Upregulated", "fdr_Downregulated", "fdr_Upregulated")
medians.out.df$Feature <- sourceMolTable[match(medians.out.df$Property, sourceMolTable$Property), "Property_formatted_short"]
medians.out.df <- medians.out.df[,c("Feature", "median_Unresponsive", "median_Downregulated", "median_Upregulated", "fdr_Downregulated", "fdr_Upregulated")]
write.table(x = medians.out.df, file = outTable, quote = FALSE, row.names = FALSE, sep = "\t")

## Compare both data frames
unchanged <- (pvals.adj.df[,2:ncol(pvals.adj.df)] > 0.05) | (medians.comp.df[,2:ncol(medians.comp.df)] == "Equal") # TRUE if non-significant or equal medians
higher <- (pvals.adj.df[,2:ncol(pvals.adj.df)] <= 0.05) & (medians.comp.df[,2:ncol(medians.comp.df)] == "Higher") # TRUE if significant and higher median
properties.heatmap.df <- data.frame(ifelse(unchanged == T, yes = 0, no = ifelse(higher == T, yes = 2, no =  1)), stringsAsFactors = F, row.names = properties)

# Colors for the heatmap
colorsH <- c("grey76", "seagreen", wes_palette("Darjeeling1")[3])
# Order the data frame to display the columns in the desired order
properties.heatmap.df <- properties.heatmap.df[,c("upregulated", "downregulated")]

# Add molecules information
properties.heatmap.df$Molecule = sourceMolTable[match(row.names(properties.heatmap.df), sourceMolTable$Property), "Molecule"]
# Order by molecule and then alphabetically
properties.hm.split = split(properties.heatmap.df, properties.heatmap.df$Molecule)
properties.hm.split = properties.hm.split[c("DNA", "RNA", "Protein")]
properties.heatmap.df = Reduce(rbind, lapply(properties.hm.split, function(x) x[order(row.names(x)),1:2]))

# Annotation colors
ann_colors = list(Molecule = c(DNA = "dodgerblue4", RNA = "deepskyblue", Protein = "#B40F20"))
# Formatted names for the properties
propsFormatted <- sourceMolTable[match(row.names(properties.heatmap.df), sourceMolTable$Property),"Property_formatted"]

# Add rownames to annotation dataframe required by pheatmap to match rows and row annotations
row.names(sourceMolTable) = sourceMolTable$Property

pheatmap(properties.heatmap.df[,1:2], 
         color = colorsH, 
         cluster_cols = F, 
         cluster_rows = F,
         border_color = "black",
         angle_col = 45, 
         labels_row = propsFormatted, 
         labels_col = c("Upregulated", "Downregulated"), 
         legend = FALSE,
         annotation_colors = ann_colors, 
         annotation_row = sourceMolTable[,colnames(sourceMolTable) == "Molecule", drop = F], 
         filename = paste0("figures/heatmapAnnotated_salt_full_properties", outputPrefix, ".pdf"), 
         cellwidth = 9,
         cellheight = 9,
         fontsize = 9,
         width = outWidth, 
         height = outHeight)




}