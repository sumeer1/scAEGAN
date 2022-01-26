# scAEGAN- Unification of Single-Cell Genomics Data by Adversarial Learning of Latent Space Correspondences 
This repository contains the online data and scAEGAN code to analyze and visualize multi-omics integration analysis, and it's downstream analysis outputs. Metrics are also available for quantifying outputs quality.

# scAEGAN Workflow
<img width="873" alt="scAEGAN" src="https://user-images.githubusercontent.com/70262340/150944062-c9c72e62-ee8b-41f2-8d97-8d7e8711529a.PNG">



# Required libraries for scAEGAN
```
The required libraries are inlcuded in environment.yml file. In order to install these libraries, follow the following steps:

1. Creating the conda environment. 
      - This will create and install the libraries included in the environment.yml file for training the scAEGAN.
        conda env create --prefix ./env --file environment.yml --force

2. The second step is to activate the conda envirnoment. 
       conda activate ./env
       

```
# Required libraries for evaluation

* To evaluate with the given scripts in the Analysis folder on the scAEGAN generated data following libraries are required: 
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
1. Training the autoencoder with the given parameters to get the latent representation by running
```
python AE.py --input_file --output_file  --batch_size  --epochs --dropout --learning_rate
```
2.Training the cyclegan with the given parameters on latent representations obtained from the Autoencoder by running

```
python cGANtrain.py --data_path --train_file --test_file --save_path --input_shape --batch_size --epochs
```

