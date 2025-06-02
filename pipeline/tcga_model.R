# optimise parameters
library(caret)  # For grid search
library(mltools)
library(MLmetrics)
library(pbapply)
library(randomForestSRC)
library(dplyr)

setwd("~/Desktop/AML_datathon/")

load(file = "model.RData")
load(file = "./results/feature_results/features.RData")

#transcriptional-signatures
sample_acts_pathway<-readRDS("data/tcga/patient_pathway_activities_TCGA.rds")
sample_acts_tf<-readRDS("data/tcga/patient_tf_activities_TCGA.rds")
#mutations
mutation_matrix <- readRDS("data/tcga/TCGA_gene_mutation_matrix_subset.RDS")
#Cibersort
cibersort_matrix<-readRDS("data/tcga/Cibersort_TCGA.rds")
#Karyotype
karyo<-read.delim("data/tcga/TCGA_ELN2022_classification_feature_matrix.tsv",row.names = "ELN2022_feature")

#####process#####

#pathway activity
sample_pathway_molten<-as.data.frame(t(sample_acts_pathway))
sample_pathway_molten$sample_id<-rownames(sample_pathway_molten)
sample_pathway_molten <- reshape2::melt(sample_pathway_molten, id.vars = "sample_id")

#tf activity
sample_pathway_molten<-as.data.frame(t(sample_acts_pathway))
sample_pathway_molten$sample_id<-rownames(sample_pathway_molten)
sample_pathway_molten <- reshape2::melt(sample_pathway_molten, id.vars = "sample_id")

#mutations
mutation_matrix_df <- data.frame(sample_id = rownames(mutation_matrix), mutation_matrix)
molten_mutation_matrix <- reshape2::melt(mutation_matrix_df, id.vars = "sample_id")
molten_mutation_matrix$variable<-paste0(molten_mutation_matrix$variable, ";mut")
molten_mutation_matrix<-molten_mutation_matrix[!is.na(molten_mutation_matrix$sample_id),]

#Cibersort results
cibersort_matrix<-cibersort_matrix[,-c(16:18)]
cibersort_matrix$P.value <- NULL
cibersort_molten_matrix<-cibersort_matrix
cibersort_molten_matrix$sample_id <- rownames(cibersort_molten_matrix)
cibersort_molten_matrix <- reshape2::melt(cibersort_molten_matrix, id.vars = "sample_id")

#karyotype
karyo<-t(karyo)
karyo <- data.frame(sample_id = rownames(karyo), karyo)
molten_karyo_matrix <- reshape2::melt(karyo, id.vars = "sample_id")
molten_karyo_matrix<-molten_karyo_matrix[!is.na(molten_karyo_matrix$sample_id),]


signature_info <- rbind(data.frame(statistic="NULL", source=sample_pathway_molten$variable,
                                   condition=sample_pathway_molten$sample_id, score=sample_pathway_molten$value,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic="NULL", source=paste0('cibersort_', cibersort_molten_matrix$variable),
                                   condition=cibersort_molten_matrix$sample_id, score=cibersort_molten_matrix$value,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic=sample_acts_tf$statistic, source=paste0('TF_', sample_acts_tf$source),
                                   condition=sample_acts_tf$condition, score=sample_acts_tf$score,
                                   p_value="NULL", feature_type = "transcriptional-signatures"),
                        data.frame(statistic="NULL", source=molten_mutation_matrix$variable,
                                   condition=molten_mutation_matrix$sample_id, score=molten_mutation_matrix$value,
                                   p_value="NULL", feature_type = "mutations"),
                        data.frame(statistic="NULL", source=paste0('karyotype_', molten_karyo_matrix$variable),
                                   condition=molten_karyo_matrix$sample_id, score=molten_karyo_matrix$value,
                                   p_value="NULL", feature_type = "karyo")
)

signature_info$condition <- gsub(signature_info$condition, pattern = '\\.', replacement='-')
clean_ids <- sub("^(([^-]+-){2}[^-]+)-.*$", "\\1", signature_info$condition)
signature_info$condition <- clean_ids

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

prep_data<-function(full_feature_data, feature_set){
  colnames(full_feature_data)<-gsub("-", "_", colnames(full_feature_data))
  colnames(full_feature_data)<-gsub(";", ".", colnames(full_feature_data))
  feature_set<-gsub("-", "_", feature_set)
  feature_set<-gsub(";", ".", feature_set)
  
  colnames(full_feature_data)<-gsub("-", "_", colnames(full_feature_data))
  colnames(full_feature_data)<-gsub(";", ".", colnames(full_feature_data))
  
  return(full_feature_data[,c(feature_set)])
}

features_to_select<-selected.vars
train_selected_features <- prep_data(signature_transposed, features_to_select)

# # Compare training and test variables
# train_vars <- o$xvar.names
# test_vars <- colnames(train_selected_features)
# missing_in_test <- setdiff(train_vars, test_vars)
# extra_in_test <- setdiff(test_vars, train_vars)
# common_vars <- intersect(train_vars, test_vars)
# train_types <- sapply(o$xvar[common_vars], class)
# test_types <- sapply(train_selected_features[, common_vars, drop = FALSE], class)
# type_mismatches <- common_vars[train_types != test_types]
# cat("Missing in test data:\n"); print(missing_in_test)

train_selected_features$RUNX1.mut <- train_selected_features$karyotype_RUNX1
train_selected_features$TP53.mut <- train_selected_features$karyotype_TP53
dim(train_selected_features)
predictions <- predict(o, newdata = train_selected_features, proximity = T)
predictions

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

# Calculate average survival time per cluster
avg_surv_by_cluster <- data_to_plot %>%
  group_by(kmeans_cluster) %>%
  summarize(avg_surv = mean(TimeSurv, na.rm = TRUE)) %>%
  arrange(avg_surv)  # from lowest (worst) to highest (best)

# Create mapping from old cluster labels to new ordered labels
# Since you want 4 = worst and 1 = best, rank clusters accordingly:
# We'll assign ranks in descending order of avg_surv (highest = 1)
avg_surv_by_cluster <- avg_surv_by_cluster %>%
  mutate(new_label = rank(-avg_surv))  # Negative for descending order

# Create a named vector for mapping old -> new
cluster_map <- setNames(avg_surv_by_cluster$new_label, avg_surv_by_cluster$kmeans_cluster)

# Apply mapping to relabel clusters
data_to_plot$kmeans_cluster_ordered <- cluster_map[as.character(data_to_plot$kmeans_cluster)]
test_proximity_clusters<-data_to_plot$kmeans_cluster_ordered[-trainIndex]

data_to_plot$ELN2017 <- sample_ELNs
