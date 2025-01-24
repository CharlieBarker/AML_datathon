
# Ensure necessary libraries are loaded
library(survival)
library(survminer)


setwd("~/Desktop/AMLproject/")
filename <- "./data/Clin_RNA.rds"
als_clin <- readRDS(filename)
Proximity_Clusters<-read.csv("../AML_datathon/results/Proximity_Clusters.csv")
# Filter the dataset to include only the relevant groups (Favorable, Intermediate, Adverse)
data <- als_clin[als_clin$ELN2017 %in% c("Favorable", "Intermediate", "Adverse"),]
data<-merge(x = Proximity_Clusters, y = data,
      by.x = "SampleID", by.y = "dbgap_rnaseq_sample")

#data<-data[data$NPM1 %in% "positive",]
# Ensure the relevant columns are in the correct format
data$overallSurvival <- as.numeric(data$overallSurvival)  # Ensure survival time is numeric
data$vitalStatus <- ifelse(data$vitalStatus == "Dead", 1, 0)  # Convert "Dead" to 1, "Alive" to 0
data<-data[data$SampleID %in% trainSampleIDs,]
# Create a Surv object for survival time and event status
surv_obj <- Surv(data$overallSurvival, data$vitalStatus)

# Fit Kaplan-Meier survival curve overall
km_fit <- survfit(surv_obj ~ 1)


table(data.frame(eln=data$ELN2017, kmeans_proximity=data$Proximity_Cluster))

library(survival)
library(survminer)
library(cowplot)

# Kaplan-Meier fit stratified by Proximity_Cluster
km_fit_cluster <- survfit(surv_obj ~ data$Proximity_Cluster)
p1 <- ggsurvplot(km_fit_cluster,
                 data = data,
                 xlab = "Time (Days)",
                 ylab = "Survival Probability",
                 title = "Proximity Cluster",
                 conf.int = F,
                 risk.table = F,
                 pval = TRUE,
                 legend.title = "")

# Kaplan-Meier fit stratified by ELN2017
km_fit_eln <- survfit(surv_obj ~ data$ELN2017)
p2 <- ggsurvplot(km_fit_eln,
                 data = data,
                 xlab = "Time (Days)",
                 ylab = "Survival Probability",
                 title = "ELN 2017",
                 conf.int = F,
                 risk.table = F,
                 pval = TRUE,
                 legend.title = "")

# Combine the two plots side by side using cowplot
combined_plot <- plot_grid(p1$plot, p2$plot, ncol = 2, labels = c("A", "B"))

# Combine the risk tables below the plots
combined_table <- plot_grid(p1$table, p2$table, ncol = 2, labels = c("", ""), rel_heights = c(1))

pdf(file = "results/km_plots.pdf", width = 15, height = 10)
# Final layout: KM plots above and risk tables below
plot_grid(combined_plot, combined_table, ncol = 1, rel_heights = c(3, 1))
dev.off()
