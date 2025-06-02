# Load necessary libraries
library(dplyr)
library(tidyr)
library(caret)
library(pROC)
library(tidyverse)
library(survival)
library(randomForestSRC)
setwd("~/Desktop/AML_datathon/")
set.seed(123)  # Set seed for reproducibility


#####load data#####


load(file = "data/prelim_rfsrc.RData")

#clinical-readings
als_clin <- readRDS("data/Clin_RNA.rds")
clin.data<-readxl::read_excel("data/beataml_wv1to4_clinical.xlsx")

#transcriptional-signatures
sample_acts_pathway<-readRDS("results/patient_pathway_activities.rds")
sample_acts_tf<-readRDS("results/patient_tf_activities.rds")
#mutations
mutation_matrix <- readRDS("data/beatAML2_gene_mutation_matrix_subset.RDS")
#Cibersort
cibersort_matrix<-read.delim("results/CIBERSORTx_BeatAML.txt")
#Karyotype
karyo<-read.delim("data/beataml_ELN2022_karyotype_renamed.tsv",row.names = "ELN2022_feature")
cols<-c("CEBPA_Biallelic","consensusAMLFusions","ageAtDiagnosis","specificDxAtAcquisition","ELN2017",
        "ageAtSpecimenAcquisition", "%.Basophils.in.PB","%.Blasts.in.BM","%.Blasts.in.PB","%.Eosinophils.in.PB",
        "%.Immature.Granulocytes.in.PB","%.Lymphocytes.in.PB","%.Monocytes.in.PB","%.Neutrophils.in.PB",
        "%.Nucleated.RBCs.in.PB","ALT","AST","albumin","creatinine","fabBlastMorphology","hematocrit",
        "hemoglobin", "otherCytogenetics","plateletCount", "totalProtein","wbcCount","FLT3-ITD",
        "allelic_ratio","NPM1","RUNX1","ASXL1","TP53")


#####process#####

transcriptional_signature_names<-unique(molten_switch_scores$variable)

#pathway activity
sample_pathway_molten<-as.data.frame(t(sample_acts_pathway))
sample_pathway_molten$sample_id<-rownames(sample_pathway_molten)
sample_pathway_molten <- reshape2::melt(sample_pathway_molten, id.vars = "sample_id")

#mutations
mutation_matrix_df <- data.frame(sample_id = rownames(mutation_matrix), mutation_matrix)
molten_mutation_matrix <- reshape2::melt(mutation_matrix_df, id.vars = "sample_id")
molten_mutation_matrix$rna_id <- als_clin$dbgap_rnaseq_sample[match(molten_mutation_matrix$sample_id, als_clin$dbgap_dnaseq_sample)]
molten_mutation_matrix$variable<-paste0(molten_mutation_matrix$variable, ";mut")
molten_mutation_matrix<-molten_mutation_matrix[!is.na(molten_mutation_matrix$rna_id),]

#Cibersort results
cibersort_matrix<-cibersort_matrix[,-c(16:19)]
cibersort_molten_matrix<-cibersort_matrix
colnames(cibersort_molten_matrix)[1]<-"sample_id"
cibersort_molten_matrix <- reshape2::melt(cibersort_molten_matrix, id.vars = "sample_id")

#karyotype
karyo<-t(karyo)
karyo <- data.frame(sample_id = rownames(karyo), karyo)
molten_karyo_matrix <- reshape2::melt(karyo, id.vars = "sample_id")
molten_karyo_matrix$rna_id <- als_clin$dbgap_rnaseq_sample[match(molten_karyo_matrix$sample_id, als_clin$dbgap_dnaseq_sample)]
molten_karyo_matrix<-molten_karyo_matrix[!is.na(molten_karyo_matrix$rna_id),]

#prognostic markers
samples_to_see<-c("Favorable", "Adverse","Intermediate","FavorableOrIntermediate","IntermediateOrAdverse")
initial_diagnosis_als_clin<-als_clin[als_clin$diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis"),]

# do we stratify on adverse only??
# initial_diagnosis_als_clin<-initial_diagnosis_als_clin[initial_diagnosis_als_clin$ELN2017 %in% c("Adverse"),]

clin_info<-data.frame(vitalStatus  = initial_diagnosis_als_clin$vitalStatus[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
                      ELN2017 = initial_diagnosis_als_clin$ELN2017[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
                      TimeSurv = initial_diagnosis_als_clin$overallSurvival[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
                      RNAseqID = initial_diagnosis_als_clin$dbgap_rnaseq_sample[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see])
signature_info <- rbind(data.frame(statistic="NULL", source=sample_pathway_molten$variable,
                                   condition=sample_pathway_molten$sample, score=sample_pathway_molten$value,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic="NULL", source=paste0('cibersort_', cibersort_molten_matrix$variable),
                                   condition=cibersort_molten_matrix$sample, score=cibersort_molten_matrix$value,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic=sample_acts_tf$statistic, source=paste0('TF_', sample_acts_tf$source),
                                   condition=sample_acts_tf$condition, score=sample_acts_tf$score,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic="NULL", source=molten_mutation_matrix$variable,
                                   condition=molten_mutation_matrix$rna_id, score=molten_mutation_matrix$value,
                                   p_value="NULL", feature_type = "mutations"),
                        data.frame(statistic="NULL", source=paste0('karyotype_', molten_karyo_matrix$variable),
                                   condition=molten_karyo_matrix$rna_id, score=molten_karyo_matrix$value,
                                   p_value="NULL", feature_type = "karyo")
)
#signature_info$score <- as.character(signature_info$score)
#signature_info<-signature_info[signature_info$feature_type == 'transcriptional-signatures',]

# Convert only mutation scores to factor
signature_info$score[signature_info$feature_type == "mutations"] <- 
  as.factor(signature_info$score[signature_info$feature_type == "mutations"])
feature_types<-unique(signature_info$feature_type)

signature_info$score[signature_info$feature_type == "karyo"] <- 
  as.factor(signature_info$score[signature_info$feature_type == "karyo"])
feature_types<-unique(signature_info$feature_type)

signature_clean <- signature_info %>%
  dplyr::select(-statistic, -p_value, -feature_type)

# Step 2: Pivot so each 'source' becomes a row, and 'condition' (sample) becomes a column
signature_wide <- signature_clean %>%
  pivot_wider(names_from = condition, values_from = score)

# Step 3: Set 'source' as rownames and remove it as a column
signature_matrix <- signature_wide %>%
  column_to_rownames("source")

# Now, 'signature_matrix' has features (from 'source') as row names and samples (conditions) as columns# Transpose signature_wide with 'source' as column names
signature_transposed <- signature_wide %>%
  column_to_rownames("source") %>%  # Use 'source' as rownames
  t() %>%                           # Transpose: samples become rows
  as.data.frame() %>%
  rownames_to_column(var = "RNAseqID")  # Add sample IDs back as a column
# 3. Merge with the clin_info data (the number of rows should align now)
merged_data <- dplyr::right_join(signature_transposed, clin_info, by = "RNAseqID")

# 3. Prepare the data for modeling
# Convert ELN2017 to a factor for classification if necessary
merged_data$ELN2017 <- factor(merged_data$ELN2017)

# Keep only the rows with complete cases (no NAs)
merged_data$vitalStatus <- ifelse(merged_data$vitalStatus == "Dead", 1,
                                  ifelse(merged_data$vitalStatus == "Alive", 0, NA))
merged_data<-na.omit(merged_data)
merged_data$TimeSurv<-as.numeric(as.character(merged_data$TimeSurv))


sample_IDs<-merged_data$RNAseqID
sample_ELNs<-merged_data$ELN2017

merged_data$RNAseqID <- NULL

trainIndex <- createDataPartition(merged_data$ELN2017, p = 0.7, list = FALSE)

train_data <- merged_data[trainIndex,]
test_data <- merged_data[-trainIndex,]
test_SampleIDs<-sample_IDs[-trainIndex]
trainSampleIDs<-sample_IDs[trainIndex]

test_SampleELNs<-sample_ELNs[-trainIndex]
train_SampleELNs<-sample_ELNs[trainIndex]


# Feature selection -------------------------------------------------------

# RandomForestSRC feature selection ---------------------------------------
# randomforestSRC
train_data_rfsrc <- train_data %>% select(-ELN2017)
model <- rfsrc(Surv(TimeSurv,vitalStatus) ~ ., data = train_data_rfsrc,ntree = 1000, nodesize = 5, importance = "permute",  # Enables permutation-based VIMP
               samptype = "swr")
vimp_result <- vimp(model)
varsel_result <- var.select(model,conservative = "high")
selected.vars<-varsel_result[["md.obj"]][["topvars"]]


# Random feature selection with permutation -------------------------------

library(randomForest)
library(caret)

# Example dataset: train_data (with Status column)
train_data.feature<-train_data %>% dplyr::select(-c(TimeSurv))

train_data.feature$vitalStatus <- as.factor(train_data.feature$vitalStatus)  # Ensure Status is a factor
train.vs.eln<-train_data.feature %>% dplyr::select(c(vitalStatus,ELN2017))
train_data.feature<-train_data.feature %>% dplyr::select(-c(vitalStatus,ELN2017))
old.col.names.train<-colnames(train_data.feature)
colnames(train_data.feature)<-gsub("-", "_", colnames(train_data.feature))
col.name.train.lookup<-cbind(colnames(train_data.feature),old.col.names.train)
col.name.train.lookup<-as.data.frame(col.name.train.lookup)
colnames(col.name.train.lookup)<-c("new","old")
#col.name.train.lookup$new<-gsub("_", "", colnames(train_data.feature))
colnames(train_data.feature)<-col.name.train.lookup$new
col.name.train.lookup[,1]<-colnames(train_data.feature)
unlist.numeric<-function(x){return(unlist(x))}
train_data.feature<-apply(train_data.feature,2,unlist.numeric)
colnames(train_data.feature)<-gsub(";", ".", colnames(train_data.feature))
col.name.train.lookup<-as.data.frame(col.name.train.lookup)
train_data.feature<-cbind(train_data.feature,train.vs.eln)
train_data.feature<-train_data.feature %>% dplyr::select(-c(ELN2017))

#train_data.feature<-do.call(cbind,train_data.feature)
# Train Random Forest model on unpermuted data to get baseline feature importance
rf_model <- randomForest(vitalStatus ~ ., data = train_data.feature, ntree = 1000)
unpermuted_feature_importance <- rf_model$importance[, 1]  # Take the MeanDecreaseGini importance score

n_permutations <- 100

# Initialize vector to store counts of feature importance scores greater than unpermuted scores
feature_importance_counts <- rep(0, length(unpermuted_feature_importance))
names(feature_importance_counts) <- names(unpermuted_feature_importance)

# Loop to permute the Status column and train Random Forest on each permuted dataset
set.seed(123)  # For reproducibility
for (i in 1:n_permutations) {
  
  # Permute the Status column
  permuted_status <- sample(train_data.feature$vitalStatus)
  
  # Replace Status in train_data with permuted Status
  permuted_train_data <- train_data.feature
  permuted_train_data$vitalStatus <- permuted_status
  
  # Train Random Forest model on permuted data
  rf_model_permuted <- randomForest(vitalStatus ~ ., data = permuted_train_data, ntree = 1000)
  
  # Get feature importance scores from the permuted data
  permuted_feature_importance <- rf_model_permuted$importance[, 1]  # MeanDecreaseGini importance
  
  # Compare permuted feature importance scores to unpermuted scores
  for (j in 1:length(unpermuted_feature_importance)) {
    if (permuted_feature_importance[j] > unpermuted_feature_importance[j]) {
      feature_importance_counts[j] <- feature_importance_counts[j] + 1
    }
  }
}

# Calculate the frequency of how often each feature's importance in the permuted data exceeds the unpermuted importance
feature_importance_freq <- feature_importance_counts / n_permutations

# Select features with frequency less than or equal to 5% (i.e., those that are consistently important)
selected_features <- names(feature_importance_freq[feature_importance_freq <= 0.05])

#####save data#####

save(selected_features, feature_importance_freq, selected.vars, varsel_result, file = "./results/feature_results/features.RData")
save(merged_data, file = "./results/feature_results/features_data.RData")
save(trainIndex, file = "./results/feature_results/partition.RData")
