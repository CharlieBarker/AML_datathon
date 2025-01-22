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

#clinical-readings
als_clin <- readRDS("../AMLproject/data/Clin_RNA.rds")
cols<-c("CEBPA_Biallelic","consensusAMLFusions","ageAtDiagnosis","specificDxAtAcquisition","ELN2017",
        "ageAtSpecimenAcquisition", "%.Basophils.in.PB","%.Blasts.in.BM","%.Blasts.in.PB","%.Eosinophils.in.PB",
        "%.Immature.Granulocytes.in.PB","%.Lymphocytes.in.PB","%.Monocytes.in.PB","%.Neutrophils.in.PB",
        "%.Nucleated.RBCs.in.PB","ALT","AST","albumin","creatinine","fabBlastMorphology","hematocrit",
        "hemoglobin", "otherCytogenetics","plateletCount", "totalProtein","wbcCount","FLT3-ITD",
        "allelic_ratio","NPM1","RUNX1","ASXL1","TP53")

#transcriptional-signatures
sample_acts_pathway<-readRDS("./results/patient_pathway_activities.rds")
sample_acts_tf<-readRDS("./results/patient_tf_activities.rds")
switch_scores <- read.delim("./features/AML_danaher_plus_class_switch_scores.tsv")
molten_switch_scores<-reshape2::melt(switch_scores)
transcriptional_signature_names<-unique(molten_switch_scores$variable)

#mutations
mutation_matrix <- readRDS("./features/beatAML2_gene_mutation_matrix_subset.RDS")
mutation_matrix_df <- data.frame(sample_id = rownames(mutation_matrix), mutation_matrix)
molten_mutation_matrix <- reshape2::melt(mutation_matrix_df, id.vars = "sample_id")
molten_mutation_matrix$rna_id <- als_clin$dbgap_rnaseq_sample[match(molten_mutation_matrix$sample_id, als_clin$dbgap_dnaseq_sample)]
molten_mutation_matrix$variable<-paste0(molten_mutation_matrix$variable, ";mut")
molten_mutation_matrix<-molten_mutation_matrix[!is.na(molten_mutation_matrix$rna_id),]

#prognostic markers

samples_to_see<-c("Favorable", "Adverse")
initial_diagnosis_als_clin<-als_clin[als_clin$diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis"),]

# do we stratify on adverse only??
# initial_diagnosis_als_clin<-initial_diagnosis_als_clin[initial_diagnosis_als_clin$ELN2017 %in% c("Adverse"),]

clin_info<-data.frame(ELN2017  = initial_diagnosis_als_clin$ELN2017[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
                      RNAseqID = initial_diagnosis_als_clin$dbgap_rnaseq_sample[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see])
signature_info <- rbind(data.frame(sample_acts_pathway, feature_type = "transcriptional-signatures"),
                        data.frame(sample_acts_tf, feature_type = "transcriptional-signatures"),
                        data.frame(statistic="NULL", source=molten_switch_scores$variable,
                                   condition=molten_switch_scores$sample, score=molten_switch_scores$value,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic="NULL", source=molten_mutation_matrix$variable,
                                   condition=molten_mutation_matrix$rna_id, score=molten_mutation_matrix$value,
                                   p_value="NULL", feature_type = "mutations")
                        )
feature_types<-unique(signature_info$feature_type)

for (feature_t in feature_types) {
  cat("processing", feature_t, " \n")
  # 1. Spread the signature data so each condition becomes a feature
  signature_clean <- signature_info %>%
    filter(feature_type == feature_t) %>%   # Filter by feature_t
    dplyr::select(-statistic, -p_value, -feature_type)
  signature_wide <- signature_clean %>%
    pivot_wider(names_from = condition, values_from = score)  # Spread by condition to create features
  # 2. Transpose the signature_wide data, remove the 'source' column and set 'source' as column names in one step
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
  model_rfe <- rfe(ELN2017 ~ ., data = train_data, sizes = c(1:10), rfeControl = rfeControl(functions = rfFuncs, method = "cv", number = 30))
  # 7. Evaluate the models on the test set
  # Make predictions with the models
  pred_rfe <- predict(model_rfe, newdata = test_data)

  # Evaluate performance of the models
  # confusionMatrix(pred_rf, test_data$ELN2017)
  confusionMatrix(pred_rfe$pred, test_data$ELN2017)

  # Save all necessary data for model
  saveRDS(train_data, file = paste0("~/Desktop/AML_datathon/training/", feature_t, "/train_data.rds"))
  saveRDS(test_data, file = paste0("~/Desktop/AML_datathon/training/", feature_t, "/test_data.rds"))
  saveRDS(colnames(merged_data), file = paste0("~/Desktop/AML_datathon/training/", feature_t, "/feature_names.rds"))
  saveRDS(model_rfe, file = paste0("~/Desktop/AML_datathon/training/", feature_t, "/model_rfe.rds"))
}
