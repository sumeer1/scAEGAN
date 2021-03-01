# scAEGAN- Unification of Single-Cell Genomics Data by Adversarial Learning of Latent Space Correspondences 
This repository contains the online data and scAEGAN code to analyze and visualize multi-omics integration analysis, and it's downstream analysis outputs. Metrics are also available for quantifying outputs quality.

# scAEGAN Workflow
[scAEGAN.pdf](https://github.com/Sumeer1/scAEGAN/files/6061815/scAEGAN.pdf)



# Required libraries for scAEGAN
```
tensorflow==1.14.0
keras==2.1.0

```
# Required libraries for evaluation
```
Seurat
clusteval
scclusteval
pheatmap
ggplot2
cowplot
pdfCluster
cluster
umap
Rtsne

```
# Basic Usage
Training the cyclegan on latent representations can be perfomed by running

```
python cGANtrain.py
```

with given training parameters
