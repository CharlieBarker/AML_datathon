# Load necessary libraries
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(ggalluvial)
library(wesanderson)
setwd("~/Desktop/AMLproject/")
filename <- "./data/Clin_RNA.rds"
als_clin <- readRDS(filename)
Proximity_Clusters<-read.csv("../AML_datathon/results/Proximity_Clusters.csv")
# Filter the dataset to include only the relevant groups (Favorable, Intermediate, Adverse)
data <- als_clin[als_clin$ELN2017 %in% c("Favorable", "Intermediate", "Adverse"),]
data<-merge(x = Proximity_Clusters, y = data,
            by.x = "SampleID", by.y = "dbgap_rnaseq_sample")

data#<-data[data$SampleID %in% test_SampleIDs,]


# Function to perform Fisher's exact test for overlap
perform_fisher_test <- function(data, variable1, variable2, dependent_variable) {
  if (length(unique(data[[dependent_variable]]))<=1) {
    return("NULL")
  }
  # Create contingency table
  contingency_table <- table(data[[variable1]], data[[dependent_variable]])
  # If expected frequencies are too low, use Fisher's exact test
  fisher_result <- fisher.test(contingency_table, simulate.p.value = T)
  return(fisher_result$p.value)
}

# Define independent variables
independent_var1 <- "Proximity_Cluster"
independent_var2 <- "Proximity_Cluster"

# Get dependent variable names (excluding independent variables)
dependent_vars <- c("cohort", "consensus_sex", "inferred_sex", "reportedRace", "reportedEthnicity", "inferred_ethnicity",
                    "centerID", "consensusAMLFusions", "isRelapse", "isDenovo", "isTransformed",
                    "priorMalignancyNonMyeloid", "priorMalignancyType", "specimenType", "consensunAMLFusions",
                    "FLT3-ITD", "NPM1")

# Initialize results data frame
results <- data.frame(dependent_variable = dependent_vars, p_value = NA)

# Loop through dependent variables and perform Fisher's test
for (i in 1:length(dependent_vars)) {
  p_val <- perform_fisher_test(data, independent_var1, independent_var2, dependent_vars[i])
  results$p_value[i] <- p_val
}

# Adjust p-values for multiple testing
results$adjusted_p_value <- p.adjust(results$p_value, method = "BH") # Benjamini-Hochberg correction

# Filter for significant results
significant_results <- results %>%
  filter(!is.na(adjusted_p_value), adjusted_p_value < 0.05) %>%
  dplyr::arrange(adjusted_p_value)

# Print significant results
print(significant_results)


alluvial_data<-data.frame(hack.AML.score=paste0("hackAML", data$Proximity_Cluster),
                          ELN.score=data$ELN2017)

alluvial_data <- alluvial_data %>%
  group_by(hack.AML.score, ELN.score) %>%
  dplyr::summarise(freq = n(), .groups = "drop")

alluvial_data$hack.AML.score <- factor(alluvial_data$hack.AML.score,
                                       levels = c("hackAML2", "hackAML1", "hackAML3", "hackAML4"))
alluvial_data$ELN.score <- factor(alluvial_data$ELN.score,
                                       levels = c("Favorable", "Intermediate", "Adverse"))

library(survival)
library(survminer)
library(dplyr)
data$overallSurvival <- as.numeric(data$overallSurvival)  # Ensure survival time is numeric
data$vitalStatus <- ifelse(data$vitalStatus == "Dead", 1, 0)  #
# Ensure factors
data$hackAML <- factor(as.character(data$Proximity_Cluster),
                       levels = c("2", "1", "3", "4"))
data$ELN.score <- factor(data$ELN2017,
                         levels = c("Favorable", "Intermediate", "Adverse"))
data$shared <- paste0(data$hackAML, "__", data$ELN.score)

# Relevel the "shared" variable to set "2__Favorable" as the reference category
data$shared <- relevel(factor(data$shared), ref = "2__Favorable")

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
data_for_pie<-data.frame(shared = paste0(data$Proximity_Cluster, "__", data$ELN.score),
                         `FLT3-ITD`=data$`FLT3-ITD`,
                         NPM1=data$NPM1)

library(ggplot2)
library(dplyr)
library(tidyr)


# Count occurrences of FLT3.ITD for each 'shared' value
data_for_pie$flt3_npm1 <- paste0(data_for_pie$FLT3.ITD, data_for_pie$NPM1)
data_counts <- data_for_pie %>%
  group_by(shared, flt3_npm1) %>%
  summarise(count = n(), .groups = "drop")


pal <- wes_palette("Zissou1", 100, type = "continuous")

ggplot(as.data.frame(combined_data),
       aes(y = freq, axis1 = ELN.score, axis2 = hack.AML.score)) +
  geom_alluvium(aes(fill = statistic), width = 1/12) +
  theme_nothing() +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 10) +  # Adjust size here
  scale_fill_gradientn(colours = pal, name = "Hazard Statistic Value") +  # Adding a title for the fill legend
  ggtitle("Comparison of ELN to hackAML stratifications") +
  theme(legend.position = "right",
        legend.title = element_text(size = 16),  # Adjust the legend title size
        legend.text = element_text(size = 14),   # Adjust the legend text size
        legend.key.size = unit(1.5, "cm"),       # Adjust size of the legend boxes
        plot.margin = margin(0, 0, 0, 0))       # Remove padding around the plot

# Create pie chart for each 'shared' value
ggplot(data_counts, aes(x = "", y = count, fill = flt3_npm1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Create pie chart
  facet_wrap(~ shared) +  # Create a separate pie chart for each 'shared' value
  theme_void() +  # Remove axis and background grid
  ggtitle("Distribution of FLT3.ITD across shared values")  # Add title

dev.off()
