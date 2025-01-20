# Load required libraries efficiently
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(org.Hs.eg.db)

# Set working directory
setwd("~/Desktop/AML_datathon//")

# Load data (keep for real analysis or remove subset if needed)
filename <- "../AMLproject/data/Norm_RNA_counts.rds"
als_clin <- readRDS(filename)

# Optional: subset to a smaller portion for testing (remove for real analysis)
als_clin <- als_clin#[, 1:10]

# Get gene names from ENSEMBL to SYMBOL
gene_names <- mapIds(org.Hs.eg.db, keys = rownames(als_clin), column = "SYMBOL", keytype = "ENSEMBL")

# Make gene names unique and clean the dataset
gene_names[is.na(gene_names)] <- rownames(als_clin)[is.na(gene_names)]
als_clin_clean <- als_clin[!is.na(gene_names), ]
rownames(als_clin_clean) <- make.unique(gene_names)

# Create results directory if not exists
results_dir <- "./results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Retrieve pathway and transcription factor networks
pathway_net <- decoupleR::get_progeny(organism = 'human', top = 500)
tf_net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)

# Run pathway network analysis (MLM)
sample_acts_pathway <- decoupleR::run_mlm(mat = als_clin_clean,
                                          net = pathway_net,
                                          .source = 'source',
                                          .target = 'target',
                                          .mor = 'weight',
                                          minsize = 5)

# # Save pathway analysis results
saveRDS(sample_acts_pathway, file.path(results_dir, "patient_pathway_activities.rds"))

# Run transcription factor network analysis (ULM)
sample_acts_tf <- decoupleR::run_ulm(mat = als_clin_clean,
                                     net = tf_net,
                                     .source = 'source',
                                     .target = 'target',
                                     .mor = 'mor',
                                     minsize = 5)

# # Save transcription factor analysis results
saveRDS(sample_acts_tf, file.path(results_dir, "patient_tf_activities.rds"))
