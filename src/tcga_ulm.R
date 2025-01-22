# Load required libraries efficiently
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# Set working directory
setwd("/SAN/colcc/datathon25/teams/team8/AML/git-scripts")

# Load data (keep for real analysis or remove subset if needed)
filename <- "../../TCGA_AML_cbioportal/laml_tcga_pub/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
als_clin <- read.delim(filename)
als_clin<-als_clin[als_clin$Hugo_Symbol != "",]
als_clin$Entrez_Gene_Id <- NULL
# Optional: subset to a smaller portion for testing (remove for real analysis)
als_clin <- als_clin[, 1:10]
rownames(als_clin) <- make.unique(als_clin$Hugo_Symbol)
als_clin$Hugo_Symbol <- NULL
als_clin<-als_clin[complete.cases(als_clin),]
# Create results directory if not exists
results_dir <- "./TCGA_test"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Retrieve pathway and transcription factor networks
pathway_net <- decoupleR::get_progeny(organism = 'human', top = 500)
tf_net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)

# Run pathway network analysis (MLM)
sample_acts_pathway <- decoupleR::run_mlm(mat = als_clin,
                                          net = pathway_net,
                                          .source = 'source',
                                          .target = 'target',
                                          .mor = 'weight',
                                          minsize = 5)

# # Save pathway analysis results
saveRDS(sample_acts_pathway, file.path(results_dir, "patient_pathway_activities.rds"))

# Run transcription factor network analysis (ULM)
sample_acts_tf <- decoupleR::run_ulm(mat = als_clin,
                                     net = tf_net,
                                     .source = 'source',
                                     .target = 'target',
                                     .mor = 'mor',
                                     minsize = 5)

# # Save transcription factor analysis results
saveRDS(sample_acts_tf, file.path(results_dir, "patient_tf_activities.rds"))
