###################
#
# 
# Multinomial logistic regression and random forest modelling of human stress genes in terms of genes features
#
#
###################

# Loading libraries
library(nnet)
library(randomForest)


conds <- "salt_HeLa:unresponsive,upregulated,downregulated"
# Read data
stress.properties.df <- read.delim("results/human/humanStressConsensus_Properties_Full_Table.tsv", stringsAsFactors = F)

properties <- colnames(stress.properties.df)[4:ncol(stress.properties.df)]

# Count number of complete observations for each property
nObservations <- plyr::ldply(colSums(!sapply(stress.properties.df[, 4:ncol(stress.properties.df)], is.na)))
colnames(nObservations) <- c("Property", "Observations")
prop.names.format = read.delim("data/original_data/features_info.txt", stringsAsFactors = F)
nObservations$Property <- prop.names.format[match(nObservations$Property, prop.names.format$Property), "Property_formatted"]

# Write table 
write.table(x = nObservations, 
            file = "results/human/numberOfObservations_per_feature.tab", 
            sep = "\t", 
            quote = F, 
            row.names = F)


# In principle, we are not filtering any variable so here I am assigning stress.properties.df w/o any change
stress.properties.df.filtered <- stress.properties.df
# Same for properties
propToInclude <- properties
# Filtering properties data frame to have the small subset used in humans
prop.names.format.filtered <- prop.names.format[prop.names.format$Property %in% propToInclude,]


# Scale data
# Omit NA to have only complete observations
# Not using NA omit here bc some stress genes categories could be NA meaning that they will be filtered wrongly
stress.properties.df.filtered.nar <- stress.properties.df.filtered[which(rowSums(is.na(stress.properties.df.filtered[,4:ncol(stress.properties.df.filtered)])) == 0),]
stress.properties.df.filtered.scaled <- cbind(stress.properties.df.filtered.nar[,1:3], scale(stress.properties.df.filtered.nar[,4:ncol(stress.properties.df.filtered.nar)]))

# Modelling and representation of the output
cond <- unlist(strsplit(conds, ":"))[1]
print(cond)
depOverlap <- unlist(strsplit(unlist(strsplit(conds, ":"))[2], ","))

for (j in c(0, 100)){
  if(j == 0){
    raw.data <- stress.properties.df.filtered.nar
  } else {
    set.seed(12345)
    newIds <- c(sample(x = which(stress.properties.df.filtered.nar[,cond] == "unresponsive"), size = j), which(stress.properties.df.filtered.nar[,cond] != "unresponsive"))
    raw.data <- stress.properties.df.filtered.nar[newIds, ]
  }
  # Multinomial logistic regression
  # Raw data - Leave-one out - AIC
  # Remove genes with NA (16 genes)
  raw.data[,paste0(cond, "_2")] <- factor(raw.data[,cond])
  raw.data[,paste0(cond,"_2")] <- relevel(raw.data[,paste0(cond,"_2")], ref = "unresponsive")
  raw.data = raw.data[!is.na(raw.data[,paste0(cond, "_2")]),]

  # Formula for the multivariate approach
  raw.f <- as.formula(paste0(cond, "_2", " ~ ", paste0(propToInclude, collapse = " + ")))
  # Multinomial test
  raw.test <- multinom(raw.f, data = raw.data)
  raw.test.res <- summary(raw.test)
  # Compute z-score and p-value and write output to a table
  z <- raw.test.res$coefficients/raw.test.res$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  raw.test.out <- as.data.frame(rbind(raw.test.res$coefficients, raw.test.res$standard.errors, p))
  row.names(raw.test.out)[grep(".1", row.names(raw.test.out))] <- gsub(".1", "_standardErrors", row.names(raw.test.out)[grep(".1", row.names(raw.test.out))], fixed = T)
  row.names(raw.test.out)[grep(".2", row.names(raw.test.out))] <- gsub(".2", "_pValue", row.names(raw.test.out)[grep(".2", row.names(raw.test.out))], fixed = T)
  write.table(x = raw.test.out,
              file = paste0("results/human/multinomialTestSummary_raw_", cond, "_random_",j, ".tab"), 
              quote = F,
              sep = "\t",
              col.names = NA)
  
  # AIC
  raw.AIC.full <- raw.test$AIC
  raw.AIC.leave.one <- lapply(propToInclude, function(x){
    raw.f.leave.one <- paste0(cond, "_2", " ~ ", paste0(grep(x, propToInclude, invert = T, value = T), collapse = " + "))
    raw.test.leave.one <- multinom(raw.f.leave.one, data = raw.data)
    raw.test.leave.one$AIC
  })
  names(raw.AIC.leave.one) <- propToInclude
  raw.AIC.leave.one$full <- raw.AIC.full
  raw.AIC.df <- plyr::ldply(raw.AIC.leave.one, .id = "Property")
  colnames(raw.AIC.df)[2] = "AIC"
  # Correct manually variables for the correct format
  raw.AIC.df$Property_Formatted <- paste0(" ", prop.names.format[match(raw.AIC.df$Property, prop.names.format$Property), "Property_formatted_short"], " ")
  raw.AIC.df[raw.AIC.df$Property == "full", "Property_Formatted"] <- " Full Model "
  # Order levels
  raw.AIC.df$Property_Formatted <- factor(raw.AIC.df$Property_Formatted, levels = rev(raw.AIC.df$Property_Formatted))
  # Save data frame from latter plotting
  write.table(raw.AIC.df, paste0("results/human/dataAIC_", cond, "_random_", j,".tsv"), sep = "\t", quote = F, row.names = F)
  
  
  # Scaled data - Coefficient of regression comparison
  if(j == 0){
    scaled.data <- stress.properties.df.filtered.scaled
  } else {
    set.seed(12345)
    newIds <- c(sample(x = which(stress.properties.df.filtered.scaled[,cond] == "unresponsive"), size = j), which(stress.properties.df.filtered.scaled[,cond] != "unresponsive"))
    scaled.data <- stress.properties.df.filtered.scaled[newIds, ]
  }
  # Remove genes with NA
  scaled.data[,paste0(cond, "_2")] <- factor(scaled.data[,cond])
  scaled.data[,paste0(cond,"_2")] <- relevel(scaled.data[,paste0(cond,"_2")], ref = "unresponsive")
  scaled.data = scaled.data[!is.na(scaled.data[,paste0(cond, "_2")]),]
  scaled.f <- as.formula(paste0(cond, "_2", " ~ ", paste0(propToInclude, collapse = " + ")))
  # Multinomial test
  scaled.test <- multinom(scaled.f, data = scaled.data)
  scaled.test.res <- summary(scaled.test)
  # Compute z-score and p-value and write output to a table
  z <- scaled.test.res$coefficients/scaled.test.res$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  scaled.test.out <- as.data.frame(rbind(scaled.test.res$coefficients, scaled.test.res$standard.errors, p))
  row.names(scaled.test.out)[grep(".1", row.names(scaled.test.out))] <- gsub(".1", "_standardErrors", row.names(scaled.test.out)[grep(".1", row.names(scaled.test.out))], fixed = T)
  row.names(scaled.test.out)[grep(".2", row.names(scaled.test.out))] <- gsub(".2", "_pValue", row.names(scaled.test.out)[grep(".2", row.names(scaled.test.out))], fixed = T)
  write.table(x = scaled.test.out,
              file = paste0("results/human/multinomialTestSummary_scaled_Properties_", cond,"_random_", j,".tab"), 
              quote = F,
              sep = "\t",
              col.names = NA)
  
  # Accuracy of multinomial logistic regression: Split the data into training and validation set
  # Training Set : Validation Set = 70 : 30 (random)
  # Test
  set.seed(12345)
  train <- sample(nrow(raw.data), 0.7*nrow(raw.data), replace = F)
  TrainSet <- raw.data[train,]
  ValidSet <- raw.data[-train,]
  # Create a Multinomial Logistic Regression model with default parameters
  model1 <- multinom(raw.f, data = TrainSet)
  # Predicting on Validation set
  predValid <- predict(model1, ValidSet, type = "class")
  # Checking classification accuracy
  mean(predValid == ValidSet[, paste0(cond, "_2")])         
  # ROC curve assesment
  # Calculate the probability of new observations belonging to each class
  # prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
  prediction_for_roc_curve <- predict(model1, ValidSet, type ="prob")
  # Save for plotting
  save(prediction_for_roc_curve, ValidSet, file = paste0("data/rdata/roc_curve_multinomial_human_random_",j,".rda"))
  
  ## Univariate modelling
  # For each property, do a model and store the output (AIC)
  univariate.measures = lapply(propToInclude, function(ppty){
    print(ppty)
    f <- as.formula(paste0(cond, "_2", " ~ ", ppty))
    raw.test <- multinom(f, data = raw.data)
    raw.aic = raw.test$AIC
    return(raw.aic)
  })
  
  names(univariate.measures) <- propToInclude
  # AIC  
  univariate.aic = plyr::ldply(univariate.measures)
  colnames(univariate.aic) = c("Property", "AIC")
  univariate.aic$Property_Formatted <- paste0(" ", prop.names.format[match(univariate.aic$Property, prop.names.format$Property), "Property_formatted_short"], " ")
  # Save data frame from latter plotting
  write.table(univariate.aic, paste0("results/human/dataAIC_univariate_", cond,"_random_",j,".tsv"), sep = "\t", quote = F, row.names = F)
  
  ### Random forest: Accuracy
  # Create a Random Forest model with default parameters
  set.seed(12345)
  model1 <- randomForest(raw.f, data = TrainSet, importance = T)
  # Predicting on Validation set
  predValid <- predict(model1, ValidSet, type = "class")
  # Checking classification accuracy
  mean(predValid == ValidSet[, paste0(cond, "_2")])
  # ROC curve assesment
  # Calculate the probability of new observations belonging to each class
  # prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
  prediction_for_roc_curve <- predict(model1, ValidSet, type ="prob")
  # Save for plotting
  save(prediction_for_roc_curve, ValidSet, file = paste0("data/rdata/roc_curve_randomForest_human_random_",j,".rda"))
  
  ### Random forest: Importance
  set.seed(12345)
  fullModel1 <- randomForest(raw.f, data = raw.data, importance = T)
  write.table(importance(fullModel1), file = paste0("results/human/modelImportance_randomForest", cond, "_random_",j,".tab"),
              quote = F, 
              sep = "\t", 
              row.names = T, 
              col.names = NA)  
  
}
