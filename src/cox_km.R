
# Ensure necessary libraries are loaded
library(survival)
library(survminer)


setwd("~/Desktop/AMLproject/")
filename <- "./data/Clin_RNA.rds"
als_clin <- readRDS(filename)

# Filter the dataset to include only the relevant groups (Favorable, Intermediate, Adverse)
data <- als_clin[als_clin$ELN2017 %in% c("Favorable", "Intermediate", "Adverse"),]
#data<-data[data$NPM1 %in% "positive",]
# Ensure the relevant columns are in the correct format
data$overallSurvival <- as.numeric(data$overallSurvival)  # Ensure survival time is numeric
data$vitalStatus <- ifelse(data$vitalStatus == "Dead", 1, 0)  # Convert "Dead" to 1, "Alive" to 0

# Create a Surv object for survival time and event status
surv_obj <- Surv(data$overallSurvival, data$vitalStatus)

# Fit Kaplan-Meier survival curve overall
km_fit <- survfit(surv_obj ~ 1)

# Kaplan-Meier curve stratified by ELN 2017 risk group
km_fit_stratified <- survfit(surv_obj ~ data$ELN2017)

# Plot the stratified Kaplan-Meier curve with confidence intervals and risk table
ggsurvplot(km_fit_stratified,
           data = data,         # Correctly pass the data to ggsurvplot
           xlab = "Time (Months)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by ELN 2017",
           conf.int = TRUE,      # Show confidence intervals
           risk.table = TRUE,    # Show the risk table with number at risk
           pval = TRUE,          # Display p-value for log-rank test
           legend.title = "ELN 2017 Classification")  # Explicitly set colors using palette




# Perform a log-rank test to compare survival curves between the three ELN2017 groups
surv_obj <- Surv(data$overallSurvival, data$vitalStatus)

# Fit Kaplan-Meier curves stratified by ELN2017
km_fit_stratified <- survfit(surv_obj ~ data$ELN2017)

# Perform log-rank test to compare survival across the ELN2017 groups
log_rank_test <- survdiff(surv_obj ~ data$ELN2017)
log_rank_test



# Ensure is_relapse is a binary variable (1 for relapse, 0 for no relapse)
data$isDenovo <- ifelse(data$isTransformed == T, 1, 0)  # If is_relapse is categorical, adjust as needed

# Fit a Cox proportional hazards model with interaction between ELN2017 and is_relapse
cox_model <- coxph(surv_obj ~ data$ELN2017 * data$isDenovo, data = data)

# Summary of the Cox model
summary(cox_model)
