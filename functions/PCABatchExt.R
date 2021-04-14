# Function to produce PCA plot from prepared data from prcomp
## Data preparation: prcomp(data)$x > as.data.frame > group column --> data frame
## x = prepared data frame
## groupName1 = name of the column containing group information in x (i.e Control and NaCl)
## groupName2 = name of the column containing group information in x (e.g. Batch)

library(ggsci)
library(ggplot2)

PCABatchExt <- function(x, 
                groupName1,
                groupName2){
  # Standarizing the obtained data frame
  d <- data.frame(row.names(x), 
                  PC1 = x$PC1, 
                  PC2 = x$PC2, 
                  Group1 = x[, groupName1], 
                  Group2 = x[,groupName2],
                  row.names = 1)
  # x and y are 1st and 2nd component label is the sample name and colored by group
  g <- ggplot(data = d, aes(x = PC1, y = PC2, colour = Group2, shape = Group1)) +
  theme_classic(base_size = 9, base_family = "Arial") + 
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size = 2) + 
  # Scaling axes to take into account the min and max values of PC
  scale_x_continuous(limits = c(min(d$PC1 - 2), max(d$PC1 + 2))) + 
  scale_y_continuous(limits = c(min(d$PC2 - 2), max(d$PC2 + 2))) + 
  scale_color_npg() + 
  # Eliminating legend and setting up appearence of x and y axis
  theme(legend.position = "right", 
      axis.text.x = element_text(size = 9, colour = "black"),
      axis.text.y = element_text(size = 9, colour = "black"),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.title.x = element_text(size = 9, colour = "black"),
      axis.title.y = element_text(size = 9, colour = "black"), 
      legend.text = element_text(colour = "black", size = 9))
  return(g)
}

