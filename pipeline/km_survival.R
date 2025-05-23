
###Comparison of Random Forest and ELN 'Severe' Classification for Survival Prediction at 1, 2, and 3 Years

# Load required libraries
library(randomForestSRC)
library(pROC)
library(ggplot2)
library(dplyr)

# Define time points
time_points <- c( 365)  # in days
colors_model <- c("blue", "red", "green")
colors_severe <- c("skyblue", "salmon", "darkgreen")

data_for_roc <- test_selected_features
# Predict survival
pred <- predict(o, newdata = data_for_roc)

# Initialize ROC data frame
roc_df_all <- data.frame()

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

# Loop through each time point
for (i in seq_along(time_points)) {
  time_point <- time_points[i]
  idx <- which.min(abs(pred$time.interest - time_point))
  
  # --- Model-based prediction
  surv_prob <- pred$survival[, idx]
  risk_score <- 1 - surv_prob
  roc_model <- roc(response = data_for_roc$vitalStatus, predictor = risk_score, quiet = TRUE)
  auc_model <- round(auc(roc_model), 3)
  
  roc_data_model <- data.frame(
    FPR = 1 - roc_model$specificities,
    TPR = roc_model$sensitivities,
    Model = paste0("Model", " (AUC = ", auc_model, ")"),
    Type = "Model",
    Time = time_point
  )
  
  # --- ELN classifier
  roc_severe <- roc(response = data_for_roc$vitalStatus, predictor = data_for_roc$ELN_score, quiet = TRUE)
  auc_severe <- round(auc(roc_severe), 3)
  
  roc_data_severe <- data.frame(
    FPR = 1 - roc_severe$specificities,
    TPR = roc_severe$sensitivities,
    Model = paste0("ELN", " (AUC = ", auc_severe, ")"),
    Type = "ELN",
    Time = time_point
  )
  
  # --- Proximity classifier
  roc_prox <- roc(response = data_for_roc$vitalStatus, predictor = data_for_roc$proximity_score, quiet = TRUE)
  auc_prox <- round(auc(roc_prox), 3)
  
  roc_data_prox <- data.frame(
    FPR = 1 - roc_prox$specificities,
    TPR = roc_prox$sensitivities,
    Model = paste0("Proximity", " (AUC = ", auc_prox, ")"),
    Type = "Proximity",
    Time = time_point
  )
  
  # Combine
  roc_df_all <- bind_rows(roc_df_all, roc_data_model, roc_data_severe, roc_data_prox)
}

# Define color scale
all_colors <- c(colors_model, colors_severe)
names(all_colors) <- unique(roc_df_all$Time)

pdf(file = './plots/roc_vitalStatus.pdf')
# Plot
ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Model, linetype = Type)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = "ROC Curves at Different Time Points",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Legend",
    linetype = "Classifier"
  ) +
  coord_equal()
dev.off()

library(survival)
library(survminer)

# Assuming your data frame is called `data_for_roc` and contains:
# - TimeSurv: time to event or censoring
# - vitalStatus: event indicator (1 = event/death, 0 = censored)
# - ELN_score: categorical or numeric score
# - proximity_score: categorical or numeric score

# Optionally convert scores to factors for cleaner KM curves
data_for_roc$ELN_score <- factor(data_for_roc$ELN_score)
data_for_roc$proximity_score <- factor(data_for_roc$proximity_score)

# Create the survival object
surv_obj <- Surv(time = data_for_roc$TimeSurv, event = data_for_roc$vitalStatus)

pdf(file = './plots/kaplein_meir.pdf')

# --- Kaplan-Meier for ELN
fit_eln <- survfit(surv_obj ~ ELN_score, data = data_for_roc)
ggsurvplot(fit_eln, data = data_for_roc,
           title = "Kaplan-Meier Curve by ELN Score",
           legend.title = "ELN Score",
           xlab = "Time (Days)", ylab = "Survival Probability",
           risk.table = TRUE, pval = TRUE)

# --- Kaplan-Meier for Proximity
fit_prox <- survfit(surv_obj ~ proximity_score, data = data_for_roc)
ggsurvplot(fit_prox, data = data_for_roc,
           title = "Kaplan-Meier Curve by Proximity Score",
           legend.title = "Proximity Score",
           xlab = "Time (Days)", ylab = "Survival Probability",
           risk.table = TRUE, pval = TRUE)

dev.off()
