# Load necessary libraries
library(dplyr)
library(tidyr)
library(caret)
library(pROC)
library(tidyverse)
library(survival)
library(survminer)

setwd("~/Desktop/AML_datathon//")

#load data
sample_acts_pathway<-readRDS("./results/patient_pathway_activities.rds")
sample_acts_tf<-readRDS("./results/patient_tf_activities.rds")
als_clin <- readRDS("../AMLproject//data/Clin_RNA.rds")
#only DNA and RNA please
#als_clin<-als_clin[als_clin$manuscript_dnaseq == "yes" & als_clin$manuscript_rnaseq == "yes",]
switch_scores <- read.delim("./features/AML_danaher_plus_class_switch_scores.tsv")
molten_switch_scores<-reshape2::melt(switch_scores)
mutation_matrix <- readRDS("./features/beatAML2_gene_mutation_matrix_subset.RDS")
mutation_matrix_df <- data.frame(sample_id = rownames(mutation_matrix), mutation_matrix)
molten_mutation_matrix <- reshape2::melt(mutation_matrix_df, id.vars = "sample_id")
molten_mutation_matrix$rna_id <- als_clin$dbgap_rnaseq_sample[match(molten_mutation_matrix$sample_id, als_clin$dbgap_dnaseq_sample)]
molten_mutation_matrix$variable<-paste0(molten_mutation_matrix$variable, ";mut")
molten_mutation_matrix<-molten_mutation_matrix[!is.na(molten_mutation_matrix$rna_id),]

samples_to_see<-c("Alive", "Dead")
initial_diagnosis_als_clin<-als_clin[als_clin$diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis"),]

# do we stratify on adverse only??
# initial_diagnosis_als_clin<-initial_diagnosis_als_clin[initial_diagnosis_als_clin$ELN2017 %in% c("Adverse"),]

clin_info<-data.frame(ELN2017  = initial_diagnosis_als_clin$vitalStatus[initial_diagnosis_als_clin$vitalStatus %in% samples_to_see],
                      RNAseqID = initial_diagnosis_als_clin$dbgap_rnaseq_sample[initial_diagnosis_als_clin$vitalStatus %in% samples_to_see])
signature_info <- rbind(sample_acts_pathway,
                        sample_acts_tf,
                        data.frame(statistic="NULL", source=molten_switch_scores$variable,
                                   condition=molten_switch_scores$sample, score=molten_switch_scores$value,
                                   p_value="NULL"),
                        data.frame(statistic="NULL", source=molten_mutation_matrix$variable,
                                   condition=molten_mutation_matrix$rna_id, score=molten_mutation_matrix$value,
                                   p_value="NULL")
                        )

# 1. Spread the signature data so each condition becomes a feature
signature_clean <- signature_info %>%
  dplyr::select(-statistic, -p_value)
signature_wide <- signature_clean %>%
  pivot_wider(names_from = condition, values_from = score)  # Spread by condition to create features

# 1. Transpose the signature_wide data, remove the 'source' column and set 'source' as column names in one step
signature_transposed <- signature_wide %>%
  dplyr::select(-source) %>%            # Remove the 'source' column
  t() %>%                               # Transpose the data
  as.data.frame() %>%                   # Convert to data frame
  rownames_to_column(var = "RNAseqID") %>% # Set RNAseqID as a column
  setNames(c("RNAseqID", signature_wide$source))  # Set the column names to 'source' values directly
# 3. Merge with the clin_info data (the number of rows should align now)
merged_data <- dplyr::right_join(signature_transposed, clin_info, by = "RNAseqID")

# 3. Prepare the data for modeling
# Convert ELN2017 to a factor for classification if necessary
merged_data$ELN2017 <- factor(merged_data$ELN2017)
merged_data$RNAseqID <- NULL
# Keep only the rows with complete cases (no NAs)
merged_data <- merged_data[complete.cases(merged_data), ]

# 4. Create a train-test split (80%-20%)
set.seed(123)  # Set seed for reproducibility
trainIndex <- createDataPartition(merged_data$ELN2017, p = 0.8, list = FALSE)
train_data <- merged_data[trainIndex,]
test_data <- merged_data[-trainIndex,]

# 5. Train a model and perform feature selection using different methods in caret
# Set up trainControl for cross-validation
train_control <- trainControl(method = "cv", number = 30, search = "grid")
# Run feature selection using different algorithms
# For example, using Random Forest and Recursive Feature Elimination
# model_rf <- train(ELN2017 ~ ., data = train_data, method = "rf", trControl = train_control)
model_rfe <- rfe(ELN2017 ~ ., data = train_data, sizes = c(1:10), rfeControl = rfeControl(functions = rfFuncs, method = "cv", number = 30))
# 7. Evaluate the models on the test set
# Make predictions with the models
# pred_rf <- predict(model_rf, newdata = test_data)
pred_rfe <- predict(model_rfe, newdata = test_data)

# Evaluate performance of the models
# confusionMatrix(pred_rf, test_data$ELN2017)
confusionMatrix(pred_rfe$pred, test_data$ELN2017)




# Save all necessary data for model
saveRDS(train_data, file = "~/Desktop/AML_datathon/results/train_data.rds")
saveRDS(test_data, file = "~/Desktop/AML_datathon/results/test_data.rds")
# Save the column names for feature names reference
saveRDS(colnames(merged_data), file = "~/Desktop/AML_datathon/results/feature_names.rds")


