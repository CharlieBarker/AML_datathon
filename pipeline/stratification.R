# optimise parameters
library(caret)  # For grid search
library(mltools)
library(MLmetrics)
library(pbapply)

setwd("~/Desktop/CRUK_AML_PAPER")
load(file = "./features.RData")
load(file = "./features_data.RData")

set.seed(123)  # Set seed for reproducibility

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


set.seed(123)  # Set seed for reproducibility

prep_data<-function(full_feature_data, feature_set){
  colnames(full_feature_data)<-gsub("-", "_", colnames(full_feature_data))
  colnames(full_feature_data)<-gsub(";", ".", colnames(full_feature_data))
  feature_set<-gsub("-", "_", feature_set)
  feature_set<-gsub(";", ".", feature_set)
  
  return(full_feature_data[,c(feature_set, 'TimeSurv', 'vitalStatus')])
}
tune_model <- function(ntree, mtry, nodesize) {
  rfsrc(Surv(TimeSurv,vitalStatus) ~ ., data = train_selected_features, ntree = ntree, mtry = mtry, nodesize = nodesize,importance=TRUE)
}

features_to_select<-selected.vars
train_selected_features <- prep_data(train_data, features_to_select)
test_selected_features <- prep_data(test_data, features_to_select)
all_selected_features <- prep_data(merged_data, features_to_select)

tune_grid <- expand.grid(
  ntree     = c(100, 300, 500, 1000),     # More trees generally improve stability
  mtry      = c(2, 3, 5, 10),       # sqrt(p), log2(p), or fixed small values
  nodesize  = c(1, 5, 10, 15)       # Smaller nodesize → deeper trees → more complexity
)

# Install pbapply if not already installed
if (!requireNamespace("pbapply", quietly = TRUE)) {
  install.packages("pbapply")
}


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





o <- rfsrc(Surv(TimeSurv,vitalStatus) ~ ., data = train_selected_features, ntree = nt, mtry = mt, nodesize =  ns,
           importance=TRUE, proximity = T, statistics=TRUE, membership = T)


data_to_plot <- all_selected_features
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
write.csv(x = data_to_plot, file = './Results/full_patient_info.csv')

