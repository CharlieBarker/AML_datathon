#make
dir.create("~/Desktop/aml_r_env")
.libPaths("~/Desktop/aml_r_env")
install.packages(c("decoupleR", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", "ggrepel", "org.Hs.eg.db"), lib = "~/Desktop/aml_r_env")
#load
.libPaths("~/Desktop/aml_r_env")
library(dplyr)
library(ggplot2)
#compress
