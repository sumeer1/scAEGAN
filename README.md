# scAEGAN- Unification of Single-Cell Genomics Data by Adversarial Learning of Latent Space Correspondences 
This repository contains the online data and scAEGAN code to analyze and visualize multi-omics integration analysis, and it's downstream analysis outputs. Metrics are also available for quantifying outputs quality.

# scAEGAN Workflow
<img width="1219" alt="Screen Shot 2021-03-01 at 4 14 33 PM" src="https://user-images.githubusercontent.com/70262340/109501866-445ded00-7aa9-11eb-9ca0-90b44a6bf091.png">



# Required libraries for scAEGAN
```
The required libraries are inlcuded in environment.yml file. In order to install these libraries, follow the following steps:

1. Creating the conda environment. This will install the libraries included in the environment.yml file for training the scAEGAN.
       conda env create --prefix ./env --file environment.yml --force

2. The second step is to activate the conda envirnoment. 
       conda activate ./env
       

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
There are two steps for the basic usage after activating the conda environment.
1. Training the autoencoder to get the latent representation by running
```
python AE.py
```
2.Training the cyclegan on latent representations obtained from the Autoencoder by running

```
python cGANtrain.py
```

with given training parameters
