# Load necessary libraries
library(dplyr)
library(tidyr)
library(caret)
library(pROC)
setwd("~/Desktop/AML_datathon//")

sample_acts_pathway<-readRDS("./results/patient_pathway_activities.rds")
sample_acts_tf<-readRDS("./results/patient_tf_activities.rds")
als_clin <- readRDS("../AMLproject//data/Clin_RNA.rds")

samples_to_see<-c("Alive", "Dead")
initial_diagnosis_als_clin<-als_clin[als_clin$diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis"),]
clin_info<-data.frame(ELN2017 = initial_diagnosis_als_clin$vitalStatus[initial_diagnosis_als_clin$vitalStatus %in% samples_to_see],
             RNAseqID=initial_diagnosis_als_clin$dbgap_rnaseq_sample[initial_diagnosis_als_clin$vitalStatus %in% samples_to_see])
signature_info <- rbind(sample_acts_pathway,
                        sample_acts_tf)

# 1. Spread the signature data so each condition becomes a feature
signature_clean <- signature_info %>%
  dplyr::select(-statistic, -p_value)
signature_wide <- signature_clean %>%
  pivot_wider(names_from = condition, values_from = score)  # Spread by condition to create features

# 1. Transpose the signature_wide data and remove the 'source' column
signature_transposed <- signature_wide %>%
  dplyr::select(-source) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "RNAseqID")  # Set RNAseqID as a column

# 2. Set the column names to the feature names (rows from original signature_wide)
colnames(signature_transposed)[-1] <- signature_wide$source  # Set condition names as columns

# 3. Merge with the clin_info data (the number of rows should align now)
merged_data <- dplyr::right_join(signature_transposed, clin_info, by = "RNAseqID")

# 3. Prepare the data for modeling
# Convert ELN2017 to a factor for classification if necessary
merged_data$ELN2017 <- factor(merged_data$ELN2017)
merged_data$RNAseqID <- NULL

# 4. Create a train-test split (80%-20%)
set.seed(123)  # Set seed for reproducibility
trainIndex <- createDataPartition(merged_data$ELN2017, p = 0.8, list = FALSE)
train_data <- merged_data[trainIndex,]
test_data <- merged_data[-trainIndex,]

# 5. Train a model and perform feature selection using different methods in caret
# Set up trainControl for cross-validation
train_control <- trainControl(method = "cv", number = 10, search = "grid")

# Run feature selection using different algorithms
# For example, using Random Forest and Recursive Feature Elimination
model_rf <- train(ELN2017 ~ ., data = train_data, method = "rf", trControl = train_control)
model_rfe <- rfe(ELN2017 ~ ., data = train_data, sizes = c(1:10), rfeControl = rfeControl(functions = rfFuncs, method = "cv", number = 10))

# 6. View the results for feature selection comparison
print(model_rf)
print(model_rfe)

# 7. Evaluate the models on the test set
# Make predictions with the models
pred_rf <- predict(model_rf, newdata = test_data)
pred_rfe <- predict(model_rfe, newdata = test_data)

# Evaluate performance of the models
confusionMatrix(pred_rf, test_data$ELN2017)
confusionMatrix(pred_rfe$pred, test_data$ELN2017)



# Extract feature importance from Random Forest model
feature_importance <- varImp(model_rf)$importance

# Add a ranking column
feature_importance$Rank <- rank(-feature_importance$Overall)  # Higher rank for more important features

# Convert rownames to a column for plotting
feature_importance <- feature_importance %>%
  rownames_to_column(var = "Feature") %>%
  arrange(desc(Overall))  # Arrange by importance

# Plot the informativeness vs. contribution
library(ggplot2)
pdf(file = "~/Desktop/test.pdf")
ggplot(feature_importance, aes(x = reorder(Feature, -Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_point(color = "red", size = 3) +
  theme_minimal() +
  labs(
    title = "Feature Informativeness vs Contribution to Classifier",
    x = "Feature (Ranked by Importance)",
    y = "Feature Importance Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


to_stratify<-clin_info[-trainIndex,]
to_stratify$predicted_status <- ifelse(to_stratify$ELN2017 == "Dead", "predicted:Dead", "predicted:Alive")  # Convert "Dead" to 1, "Alive" to 0


filename <- "./data/Clin_RNA.rds"
als_clin <- readRDS(filename)

# Rename RNAseqID in to_stratify to match dbgap_rnaseq_sample in als_clin
als_clin_train<-als_clin[als_clin$dbgap_rnaseq_sample %in% to_stratify$RNAseqID,]
als_clin_train$predicted_status <- to_stratify$ELN2017[match(to_stratify$RNAseqID, als_clin_train$dbgap_rnaseq_sample)]

# Ensure the relevant columns are in the correct format
als_clin_train$overallSurvival <- as.numeric(als_clin_train$overallSurvival)  # Ensure survival time is numeric
als_clin_train$vitalStatus <- ifelse(als_clin_train$vitalStatus == "Dead", 1, 0)  # Convert "Dead" to 1, "Alive" to 0

# Create a Surv object for survival time and event status
surv_obj <- Surv(als_clin_train$overallSurvival, als_clin_train$vitalStatus)

# Fit Kaplan-Meier survival curve overall
km_fit <- survfit(surv_obj ~ 1)

# Kaplan-Meier curve stratified by ELN 2017 risk group
km_fit_stratified <- survfit(surv_obj ~ als_clin_train$predicted_status)

# Plot the stratified Kaplan-Meier curve with confidence intervals and risk table
ggsurvplot(km_fit_stratified,
           data = als_clin_train,         # Correctly pass the data to ggsurvplot
           xlab = "Time (Months)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by predicted status",
           conf.int = TRUE,      # Show confidence intervals
           risk.table = TRUE,    # Show the risk table with number at risk
           pval = TRUE)  # Explicitly set colors using palette

dev.off()

