library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(caret)

# Set minimum library size
minlib = 15000

# Select unfractionated samples (whole bacterial community)
# Select samples based on Blastocystis presence/absence and scale reads to minlib library size
PS %>%
  subset_samples(Blastocystis.Illumina %in% c("Positive", "Negative")) 

# Remove OTUs from the family "Blastocystidae"
otu_to_remove <- which(tax_table(PS)[, "family"] == "f__Blastocystidae")

if (length(otu_to_remove) > 0) {
  PS <- prune_taxa(taxa_names(PS)[-otu_to_remove], PS)
}

# Print number of OTUs
ntaxa(PS)

# Set pruning scale for rare OTUs
prunescale = 0.0001

# Prune rare OTUs based on mean relative abundance determined by prunescale
tax.mean <- taxa_sums(PS)/nsamples(PS)
sites.prune <- prune_taxa(tax.mean > prunescale * minlib, PS)

# Print pruned data
sites.prune

# Create a dataframe of training data with OTUs as columns and samples as rows
predictors <- t(otu_table(sites.prune))
dim(predictors)

# Create one column for the outcome/response variable
response <- as.factor(sample_data(sites.prune)$Blastocystis.Illumina)

# Combine predictors and response into one data frame
rf.data <- data.frame(response, predictors)

# Set seed for reproducibility
set.seed(2)

# Train a random forest model
erie.classify <- randomForest(response ~ ., data = rf.data, ntree = 100)
print(erie.classify)

# List names of the random forest model elements
names(erie.classify)

# Create a dataframe with predictor names and their importance scores
imp <- importance(erie.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Sort predictors by importance
imp.sort <- imp[order(imp$MeanDecreaseGini, decreasing = TRUE), ]
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]

# Plot the importance of top 20 predictors using ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying Blastocystis colonization")

# Install and import the pROC library if not already installed
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

# Calculate prediction probabilities
pred_prob <- predict(erie.classify, rf.data, type = "prob")

# Create ROC curve object
roc_curve <- roc(rf.data$response, pred_prob[,2], levels = rev(levels(rf.data$response)))

# Plot the ROC curve
plot.roc(roc_curve, col = "blue", print.auc = TRUE, main = "ROC Curve for Blastocystis Classification")

# Install and import the pdp library if not already installed
library(pdp)

# Make predictions using the random forest model
predictions <- predict(erie.classify, rf.data)

# Create a confusion matrix
conf_matrix <- confusionMatrix(predictions, rf.data$response)

# Print the confusion matrix
print(conf_matrix)

# Extract values from the confusion matrix
tn <- conf_matrix$table[1, 1]  # True negatives
fp <- conf_matrix$table[1, 2]  # False positives
fn <- conf_matrix$table[2, 1]  # False negatives
tp <- conf_matrix$table[2, 2]  # True positives

# Calculate performance metrics
precision <- tp / (tp + fp)
recall <- tp / (tp + fn)
f1_score <- 2 * (precision * recall) / (precision + recall)
specificity <- tn / (tn + fp)

# Display the performance metrics
cat("Precision:", precision, "\n")
cat("Recall (Sensitivity):", recall, "\n")
cat("F1 Score:", f1_score, "\n")
cat("Specificity:", specificity, "\n")


