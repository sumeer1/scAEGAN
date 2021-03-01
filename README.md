# scAEGAN- Unification of Single-Cell Genomics Data by Adversarial Learning of Latent Space Correspondences 
This repository contains the online data and scAEGAN code to analyze and visualize multi-omics integration analysis, and it's downstream analysis outputs. Metrics are also available for quantifying outputs quality.

# scAEGAN Workflow
<img width="1219" alt="Screen Shot 2021-03-01 at 4 14 33 PM" src="https://user-images.githubusercontent.com/70262340/109501866-445ded00-7aa9-11eb-9ca0-90b44a6bf091.png">



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
