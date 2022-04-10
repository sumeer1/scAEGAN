# Package names
packages <- c("ggplot2", "Seurat", "tidyverse", "cowplot", "clusteval", "pdfCluster", "pheatmap", "cluster", "Rtsne", "umap", "e1071")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
