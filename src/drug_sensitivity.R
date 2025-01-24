# Load necessary libraries
library(ggplot2)
library(rstatix)
library(ggpubr)
setwd("~/Desktop/AML_datathon//")

drug_Response<-read.csv("~/Downloads/AUC_Drug_response.csv")
colnames(drug_Response)<-drug_Response[1,]
drug_Response<-drug_Response[-1,-1]
proximity_clusters<-read_csv("./results/Proximity_Clusters.csv")

filtered_drug_Response <- drug_Response %>%
  filter(drug_Response$dbgap_rnaseq_sample %in% proximity_clusters$SampleID)

filtered_drug_Response$Proximity_Cluster<-proximity_clusters$Proximity_Cluster[match(filtered_drug_Response$dbgap_rnaseq_sample,proximity_clusters$SampleID)]
filtered_drug_Response$Proximity_Cluster<-as.character(filtered_drug_Response$Proximity_Cluster)
filtered_drug_Response$auc<-as.numeric(as.character(filtered_drug_Response$auc))
unique.drugs<-unique(filtered_drug_Response$inhibitor)


to_test<-data.frame(patient = filtered_drug_Response$dbgap_subject_id,
                    drug = filtered_drug_Response$inhibitor,
                    AUC = filtered_drug_Response$auc,
                    group = filtered_drug_Response$Proximity_Cluster)
# Your dataset
data <- to_test

# Step 1: Generate synthetic group data for testing (replace this with real groups)
set.seed(42)  # For reproducibility

# Step 2: Perform Kruskal-Wallis test per drug
unique_drugs <- unique(data$drug)
results_drugs <- list()  # Store p-values

for (drug in unique_drugs) {
  subset_data <- data[data$drug == drug,]

  if (length(unique(subset_data$group)) > 1) {
    test <- kruskal.test(AUC ~ as.factor(group), data = subset_data)
    results_drugs[[drug]]<-data.frame(pval=test$p.value, stat=unname(test$statistic))
  } else {
    p_values <- NA
    statistics <- NA
    results_drugs[[drug]]<-data.frame(pval=p_values, stat=statistics)
  }
}
complete_drugs<-bind_rows(results_drugs, .id="drug_name")
complete_drugs$p.adj <- p.adjust(complete_drugs$pval)

# Step 3: Adjust p-values for FDR
sig_drugs<-complete_drugs[complete_drugs$p.adj < 0.01,]
sig_to_test<-to_test[to_test$drug %in% sig_drugs$drug_name,]




# Step 1: Perform regression for each group and drug
unique_drugs <- unique(sig_to_test$drug)
unique_groups <- unique(sig_to_test$group)

coeff_matrix <- matrix(NA, nrow = length(unique_groups), ncol = length(unique_drugs),
                       dimnames = list(paste0("Group_", unique_groups), unique_drugs))
pval_matrix <- matrix(NA, nrow = length(unique_groups), ncol = length(unique_drugs),
                      dimnames = list(paste0("Group_", unique_groups), unique_drugs))

for (group in unique_groups) {
  for (drug in unique_drugs) {
    subset_data <- data[data$drug == drug,]
    subset_data$binary_group <- ifelse(subset_data$group == group, 1, 0)  # Compare group vs others

    if (length(unique(subset_data$binary_group)) > 1) {
      model <- glm(binary_group ~ AUC, data = subset_data, family = binomial())
      coeff_matrix[paste0("Group_", group), drug] <- coef(model)["AUC"]
      pval_matrix[paste0("Group_", group), drug] <- summary(model)$coefficients["AUC", "Pr(>|z|)"]
    }
  }
}

# Step 2: Adjust p-values for FDR
adjusted_pval_matrix <- apply(pval_matrix, 2, function(col) p.adjust(col, method = "fdr"))

# Create a character matrix based on the adjusted p-values
significance_matrix <- apply(adjusted_pval_matrix, c(1, 2), function(p) {
  if (is.na(p)) {
    return("")  # If the value is NA, leave it as an empty string
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.1) {
    return("*")
  } else {
    return("")  # No significance if p-value is greater than 0.1
  }
})

library(ComplexHeatmap)
pdf("~/Desktop/test.pdf", width = 12, height = 8)
Heatmap(
  coeff_matrix,  # Transpose so rows are drugs, columns are groups
  name = "Coefficients",  # Legend title
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Blue to red gradient
  cluster_rows = TRUE,  # Cluster rows (drugs)
  cluster_columns = TRUE,  # Cluster columns (groups)
  na_col = "darkgrey",  # Color for NA values
  row_title = "Drugs",  # Label for rows
  column_title = "Groups",  # Label for columns
  rect_gp = gpar(col = "white", lwd = 2),  # Cell border settings
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Add significance labels from the significance matrix
    sig_label <- t(significance_matrix)[j, i]
    if (sig_label != "") {
      grid.text(sig_label, x, y, gp = gpar(fontsize = 8, col = "black"))
    }
  },
  heatmap_width = unit(15, "cm"),  # Adjust width
  heatmap_height = unit(10, "cm"),  # Adjust height
  column_title_gp = gpar(fontsize = 12),  # Font size for column title
  row_title_gp = gpar(fontsize = 12)  # Font size for row title
)

# Beautified boxplot with jitter
library(ggplot2)

# Subset data for Panobinostat
Panobinostat <- filtered_drug_Response[filtered_drug_Response$inhibitor %in% "Palbociclib",]

# Create the plot
p <- ggplot(Panobinostat, aes(x = Proximity_Cluster, y = auc, fill = Proximity_Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") + # Boxplot with no outliers shown
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 2, color = "darkblue") + # Jittered points
  cowplot::theme_cowplot() + # Minimal theme for a clean look
  scale_fill_brewer(palette = "Set3") + # Set3 palette for vibrant colors
  labs(title = "Palbociclib Response by Proximity Cluster",
       x = "Proximity Cluster",
       y = "AUC",
       fill = "Cluster") + # Labels for axes and legend
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Center and style title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top"
  )

# Save the plot as a PDF
print(p)

dev.off()


