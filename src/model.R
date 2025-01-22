# Load necessary libraries
library(dplyr)
library(caret)
library(pROC)
library(tidyverse)

# Load saved data from feature selection
train_data <- readRDS("~/Desktop/AML_datathon/results/train_data.rds")
test_data <- readRDS("~/Desktop/AML_datathon/results/test_data.rds")
merged_data_columns <- readRDS("~/Desktop/AML_datathon/results/feature_names.rds")

cols<-c("CEBPA_Biallelic","consensusAMLFusions","ageAtDiagnosis","specificDxAtAcquisition","ELN2017",
        "ageAtSpecimenAcquisition", "%.Basophils.in.PB","%.Blasts.in.BM","%.Blasts.in.PB","%.Eosinophils.in.PB",
        "%.Immature.Granulocytes.in.PB","%.Lymphocytes.in.PB","%.Monocytes.in.PB","%.Neutrophils.in.PB",
        "%.Nucleated.RBCs.in.PB","ALT","AST","albumin","creatinine","fabBlastMorphology","hematocrit",
        "hemoglobin", "otherCytogenetics","plateletCount", "totalProtein","wbcCount","FLT3-ITD",
        "allelic_ratio","NPM1","RUNX1","ASXL1","TP53")
# Get the selected features from RFE
# Extract feature importance from Random Forest model
feature_importance <- data.frame(varImp(model_rfe))
# Add a ranking column
feature_importance$Rank <- rank(-feature_importance$Overall)  # Higher rank for more important features
# Convert rownames to a column for plotting
feature_importance <- feature_importance %>%
  rownames_to_column(var = "Feature") %>%
  arrange(desc(Overall))  # Arrange by importance

feature_importance$Feature <- gsub("`", "", feature_importance$Feature)
selected_features <- feature_importance$Feature[1:30]
# Reduce the train and test datasets to include only selected features and target
train_data_reduced <- train_data[, c(selected_features, "ELN2017")]
test_data_reduced <- test_data[, c(selected_features, "ELN2017")]

# Define a list of models to train
model_list <- c("rf", "glm", "svmRadial")

# Set up trainControl for consistent cross-validation
train_control <- trainControl(
  method = "cv",          # Cross-validation
  number = 10,            # 10-fold CV
  classProbs = TRUE,      # Needed for AUC calculation
  summaryFunction = twoClassSummary  # Use AUC as the metric
)

# Initialize a list to store trained models and results
models <- list()
roc_curves <- list()

# Train models using the selected features
for (model in model_list) {
  set.seed(123)  # For reproducibility
  models[[model]] <- train(
    ELN2017 ~ .,
    data = train_data_reduced,
    method = model,
    trControl = train_control,
    metric = "ROC"  # Optimize based on AUC
  )
}

# Evaluate models on the test data
library(pROC)

for (model in names(models)) {
  # Get predicted probabilities for the "Dead" class
  test_probs <- predict(models[[model]], newdata = test_data_reduced, type = "prob")[, "Dead"]

  # Calculate AUC
  roc_obj <- roc(test_data_reduced$ELN2017, test_probs, levels = rev(levels(test_data_reduced$ELN2017)))
  auc_value <- auc(roc_obj)

  # Store AUC and ROC curve
  roc_curves[[model]] <- list(roc = roc_obj, auc = auc_value)
  cat(sprintf("Model: %s, AUC: %.3f\n", model, auc_value))
}

# Plot the ROC curves for all models
library(ggplot2)
pdf(file = "~/Desktop/roc_curves.pdf")

ggplot(feature_importance[feature_importance$Feature %in% selected_features,], aes(x = reorder(Feature, -Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_point(color = "red", size = 3) +
  theme_minimal() +
  labs(
    title = "Feature Informativeness vs Contribution to Classifier",
    x = "Feature (Ranked by Importance)",
    y = "Feature Importance Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot the first ROC curve (Random Forest)
plot(
  roc_curves[["rf"]]$roc,
  col = "blue",
  lwd = 2,
  main = "ROC Curves for Classifiers",
  xlab = "False Positive Rate (1 - Specificity)",
  ylab = "True Positive Rate (Sensitivity)"
)

# Add additional ROC curves
model_colors <- c("blue", "red", "green", "purple")  # Predefined colors for consistency
model_index <- 1

for (model in names(roc_curves)[-1]) {
  model_index <- model_index + 1
  plot(
    roc_curves[[model]]$roc,
    col = model_colors[model_index],
    add = TRUE,
    lwd = 2
  )
}

# Add a legend with model names and AUC values
legend_labels <- sapply(names(roc_curves), function(model) {
  paste0(model, " (AUC=", round(roc_curves[[model]]$auc, 3), ")")
})
legend(
  "bottomright",                     # Position of the legend
  legend = legend_labels,            # Labels with AUC values
  col = model_colors[1:length(roc_curves)],  # Use the corresponding colors
  lwd = 2                            # Line width
)


dev.off()


