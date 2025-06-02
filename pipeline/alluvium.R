# Load necessary libraries
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(ggalluvial)
library(wesanderson)
library(survival)
library(survminer)
library(dplyr)

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

alluvial_data<-data.frame(hack.AML.score=paste0("hackAML", data_for_roc$proximity_score),
                          ELN.score=data_for_roc$ELN_score)

alluvial_data <- alluvial_data %>%
  group_by(hack.AML.score, ELN.score) %>%
  dplyr::summarise(freq = n(), .groups = "drop")

alluvial_data$hack.AML.score <- factor(alluvial_data$hack.AML.score)
alluvial_data$ELN.score <- factor(alluvial_data$ELN.score)

data<-data_for_roc

data$overallSurvival <- as.numeric(data$TimeSurv)  # Ensure survival time is numeric

# Ensure factors
data$hackAML <- factor(as.character(data$proximity_score),
                       levels = c("2", "1", "3", "4"))
data$ELN.score <- factor(data$ELN_score)
data$shared <- paste0(data$hackAML, "__", data$ELN.score)

# Relevel the "shared" variable to set "2__Favorable" as the reference category
data$shared <- relevel(factor(data$shared), ref = "1__1")

# Create a Surv object
surv_obj <- Surv(time = data$overallSurvival, event = data$vitalStatus)

# Fit Cox proportional hazards model including both factors
cox_model <- coxph(surv_obj ~ shared, data = data)

# Extract hazard ratios and confidence intervals
hr_results <- broom::tidy(cox_model, exponentiate = TRUE)  # Converts log(HR) to HR

# Add the row for shared2__Favorable with hazard ratio of 1 and p-value of 1
new_row <- tibble(
  term = "shared2__Favorable",
  estimate = 1,
  std.error = NA,   # Not provided
  statistic = 0,   # Not provided
  p.value = 1       # p-value of 1
)

# Bind this new row to the existing hr_results
hr_results_with_control <- bind_rows(hr_results, new_row)
hr_results_with_control$term <- gsub(pattern = "shared", replacement = "hackAML", x = hr_results_with_control$term)

# Now merge with alluvial_data based on the shared variable
# First, create a new "shared" column in alluvial_data to match the levels in hr_results
alluvial_data$shared <- paste0(alluvial_data$hack.AML.score, "__", alluvial_data$ELN.score)

# Merge the two data frames on the "shared" column
combined_data <- left_join(alluvial_data, hr_results_with_control, by = c("shared" = "term"))

library(ggplot2)
library(dplyr)
library(tidyr)
pdf(file = './plots/alluvium.pdf', width = 12)

pal <- wes_palette("Zissou1", 100, type = "continuous")

ggplot(as.data.frame(combined_data),
       aes(y = freq, axis1 = ELN.score, axis2 = hack.AML.score)) +
  geom_alluvium(aes(fill = statistic), width = 1/12) +
  theme_void() +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 10) +  # Adjust size here
  scale_fill_gradientn(colours = pal, name = "Hazard Statistic Value") +  # Adding a title for the fill legend
  ggtitle("Comparison of ELN to hackAML stratifications") +
  theme(legend.position = "right",
        legend.title = element_text(size = 16),  # Adjust the legend title size
        legend.text = element_text(size = 14),   # Adjust the legend text size
        legend.key.size = unit(1.5, "cm"))       # Remove padding around the plot

dev.off()