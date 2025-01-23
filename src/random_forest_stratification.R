# Load necessary libraries
library(dplyr)
library(tidyr)
library(caret)
library(pROC)
library(tidyverse)
library(survival)
library(survminer)
library(randomForestSRC)
setwd("~/Desktop/AML_datathon//")

set.seed(123)  # Set seed for reproducibility


#load data

#clinical-readings
als_clin <- readRDS("/Users/lourdes/Documents/Postdoc_2023_Enver_lab/CRUK_AML/Data/Clin_RNA.rds")
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

samples_to_see<-c("Favorable", "Adverse","Intermediate","FavorableOrIntermediate","IntermediateOrAdverse")
initial_diagnosis_als_clin<-als_clin[als_clin$diseaseStageAtSpecimenCollection %in% c("Initial Diagnosis"),]

# do we stratify on adverse only??
# initial_diagnosis_als_clin<-initial_diagnosis_als_clin[initial_diagnosis_als_clin$ELN2017 %in% c("Adverse"),]

clin_info<-data.frame(vitalStatus  = initial_diagnosis_als_clin$vitalStatus[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
                      ELN2017 = initial_diagnosis_als_clin$ELN2017[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
                      TimeSurv = initial_diagnosis_als_clin$overallSurvival[initial_diagnosis_als_clin$ELN2017 %in% samples_to_see],
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

# 1. Spread the signature data so each condition becomes a feature
signature_clean <- signature_info %>%
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

# Keep only the rows with complete cases (no NAs)
merged_data <- merged_data[complete.cases(merged_data), ]
merged_data$vitalStatus <- ifelse(merged_data$vitalStatus == "Dead", 1,
                                  ifelse(merged_data$vitalStatus == "Alive", 0, NA))
merged_data<-na.omit(merged_data)
merged_data$TimeSurv<-as.numeric(as.character(merged_data$TimeSurv))

sample_IDs<-merged_data$RNAseqID
merged_data$RNAseqID <- NULL

trainIndex <- createDataPartition(merged_data$ELN2017, p = 0.8, list = FALSE)
train_data <- merged_data[trainIndex,]
test_data <- merged_data[-trainIndex,]
# optimise parameters
library(caret)  # For grid search
library(mltools)
library(MLmetrics)
# Define the grid of hyperparameters to search over
tune_grid <- expand.grid(
  ntree = c(20,50,100,500,1000),
  mtry = c(3,5,10,20),
  nodesize = c(3,5,10,20)
)
tune_model <- function(ntree, mtry, nodesize) {
  rfsrc(Surv(TimeSurv,vitalStatus) ~ ., data = train_data, ntree = ntree, mtry = mtry, nodesize = nodesize,importance=TRUE)
}

# Install pbapply if not already installed
if (!requireNamespace("pbapply", quietly = TRUE)) {
  install.packages("pbapply")
}

library(pbapply)

# Add progress bar to the lapply loop
results <- pblapply(1:nrow(tune_grid), function(i) {
  params <- tune_grid[i, ]
  model <- tune_model(params$ntree, params$mtry, params$nodesize)
  mean(model$err.block.rate[1])  # Collect average error rate
})
# Cross-validate over the parameter grid
# results <- lapply(1:nrow(tune_grid), function(i) {
#     params <- tune_grid[i, ]
#     model <- tune_model(params$ntree, params$mtry, params$nodesize)
#     mean(model$err.block.rate[1])  # Collect average error rate
#   })

# Select best parameter combination based on minimum error
best_params <- tune_grid[which.min(unlist(results)), ]
print(best_params)
nt<-as.numeric(best_params[1])
mt<-as.numeric(best_params[2])
ns<-as.numeric(best_params[3])

load(file = "~/Downloads/prelim_rfsrc.RData")



o <- rfsrc(Surv(TimeSurv,vitalStatus) ~ ., data = train_data, ntree = nt, mtry = mt, nodesize =  ns,
           importance=TRUE, proximity = T, statistics=TRUE, membership = T)

data_to_plot <- merged_data
# Predict on the test_data using the fitted model

predictions <- predict(o, newdata = data_to_plot, proximity = T)
# Step 1: Convert proximity matrix to dissimilarity matrix
dissimilarity_matrix <- as.dist(1 - predictions$proximity)

# Step 2: Dimensionality reduction using Multidimensional Scaling (MDS)
library(MASS)
mds <- cmdscale(dissimilarity_matrix, k = 3) # Reduce to 2 dimensions

# Step 1: Determine the optimal number of clusters using an Elbow Plot
wss <- numeric(10) # Within-cluster sum of squares
for (k in 1:10) {
  kmeans_result <- kmeans(mds, centers = k, nstart = 10)
  wss[k] <- kmeans_result$tot.withinss # Total within-cluster sum of squares
}

# Step 2: Plot the Elbow Plot
plot(1:10, wss, type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters (k)", ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Plot for Optimal Clusters")

# Step 3: Perform K-means clustering with the optimal number of clusters
optimal_k <- 4 # Choose based on the elbow plot
final_kmeans <- kmeans(mds, centers = optimal_k, nstart = 10)

# Step 4: Add cluster labels to the data
cluster_labels <- final_kmeans$cluster
data_to_plot$kmeans_cluster <- cluster_labels

to_write<-data_to_plot
write.csv(x = data.frame(Proximity_Cluster = to_write$kmeans_cluster,
                         SampleID=sample_IDs), file = "./results/Proximity_Clusters.csv")

# Ensure the relevant columns are in the correct format
train_data_plot<-data_to_plot[data_to_plot$ELN2017 %in% c("Adverse", "Favorable", "Intermediate"),]
train_data_plot$overallSurvival <- as.numeric(train_data_plot$TimeSurv)  # Ensure survival time is numeric

# Create a Surv object for survival time and event status
surv_obj <- Surv(train_data_plot$overallSurvival, train_data_plot$vitalStatus)

# Fit Kaplan-Meier survival curve overall
km_fit <- survfit(surv_obj ~ train_data_plot$kmeans_cluster)
# Plot the stratified Kaplan-Meier curve with confidence intervals and risk table
plot1<-ggsurvplot(km_fit,
                  data = train_data_plot,         # Correctly pass the data to ggsurvplot
                  xlab = "Time (Months)",
                  ylab = "Survival Probability",
                  title = "KM Survival Curve by RFsrc",
                  conf.int = TRUE,      # Show confidence intervals
                  risk.table = TRUE,    # Show the risk table with number at risk
                  pval = TRUE,          # Display p-value for log-rank test
                  legend.title = "ELN 2017 Classification")  # Explicitly set colors using palette
# Fit Kaplan-Meier survival curve overall
km_fit <- survfit(surv_obj ~ train_data_plot$ELN2017)
# Plot the stratified Kaplan-Meier curve with confidence intervals and risk table
plot2<-ggsurvplot(km_fit,
                  data = train_data_plot,         # Correctly pass the data to ggsurvplot
                  xlab = "Time (Months)",
                  ylab = "Survival Probability",
                  title = "KM Survival Curve by ELN2017",
                  conf.int = TRUE,      # Show confidence intervals
                  risk.table = TRUE,    # Show the risk table with number at risk
                  pval = TRUE,          # Display p-value for log-rank test
                  legend.title = "ELN 2017 Classification")  # Explicitly set colors using palette



# Combine plots into a single PDF
pdf(file = "~/Desktop/KM_Survival_Plots.pdf", width = 8, height = 10)
print(plot1)
print(plot2)
dev.off()

table(data.frame(eln=train_data_plot$ELN2017, kmeans_proximity=train_data_plot$kmeans_cluster))
