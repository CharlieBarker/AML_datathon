library(survival)
library(timeROC)

setwd("~/Desktop/CRUK_AML_PAPER")


load(file = "./features.RData")
load(file = "./features_data.RData")
load(file = "./partition.RData")
load(file = "./model.RData")

data_to_plot <- read.csv('./Results/full_patient_info.csv')

# Check if trainIndex exists
if (exists("trainIndex")) {
  test_proximity_clusters <- data_to_plot$kmeans_cluster_ordered[-trainIndex]
} else {
  stop("trainIndex is not defined.")
}

set.seed(123)  # Set seed for reproducibility

merged_data$RNAseqID <- NULL

test_data <- merged_data[-trainIndex,]

prep_data<-function(full_feature_data, feature_set){
  colnames(full_feature_data)<-gsub("-", "_", colnames(full_feature_data))
  colnames(full_feature_data)<-gsub(";", ".", colnames(full_feature_data))
  feature_set<-gsub("-", "_", feature_set)
  feature_set<-gsub(";", ".", feature_set)
  
  return(full_feature_data[,c(feature_set, 'TimeSurv', 'vitalStatus')])
}

features_to_select<-selected.vars
test_selected_features <- prep_data(test_data, features_to_select)

sample_IDs<-merged_data$RNAseqID
sample_ELNs<-merged_data$ELN2017
test_SampleELNs<-sample_ELNs[-trainIndex]
train_SampleELNs<-sample_ELNs[trainIndex]

data_for_roc <- test_selected_features


# Create binary predictor: Severe = 1, else = 0
# Convert ELN labels to numeric severity scores
eln_scores <- recode(
  test_SampleELNs,
  "Favorable" = 1,
  "FavorableOrIntermediate" = 1,
  "Intermediate" = 2,
  "IntermediateOrAdverse" = 2,
  "Adverse" = 3
)
data_for_roc$ELN_score <- as.numeric(eln_scores)
data_for_roc$proximity_score <- test_proximity_clusters

time_seq<-seq(0.4,.9,0.2)
# we evaluate our proximity as a prognostic biomarker.
ROC.proximity<-timeROC(T=data_for_roc$TimeSurv,
                  delta=data_for_roc$vitalStatus,
                  marker=as.numeric(data_for_roc$proximity_score),
                  cause=1,
                  weighting="marginal",
                  times=quantile(data_for_roc$TimeSurv,probs=time_seq),
                  iid=TRUE)
# we evaluate ELN as a prognostic biomarker.
ROC.ELN<-timeROC(T=data_for_roc$TimeSurv,
                  delta=data_for_roc$vitalStatus,
                  marker=as.numeric(data_for_roc$ELN_score),
                  cause=1,
                  weighting="marginal",
                  times=quantile(data_for_roc$TimeSurv,probs=time_seq),
                  iid=TRUE)

print(ROC.proximity)
print(ROC.ELN)

compare(ROC.proximity,ROC.ELN) #compute p-values of comparison tests

pdf(file = './plots/AUC_curve.pdf')
# plot AUC curve for albumin only with pointwise confidence intervals
# and simultaneous confidence bands
plotAUCcurve(ROC.proximity,conf.int=TRUE,conf.band=TRUE)
# plot AUC curve for albumin and bilirunbin with pointwise confidence intervals
plotAUCcurve(ROC.proximity,conf.int=TRUE,col="red")
plotAUCcurve(ROC.ELN,conf.int=TRUE,col="blue",add=TRUE)
legend("bottomright",c("Proxmity","ELN 2017"),col=c("red","blue"),lty=1,lwd=2)
dev.off()
