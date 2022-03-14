# scAEGAN- Unification of Single-Cell Genomics Data by Adversarial Learning of Latent Space Correspondences 
This repository contains the online data and scAEGAN code to analyze and visualize multi-omics integration analysis, and it's downstream analysis outputs. Metrics are also available for quantifying outputs quality.


# scAEGAN Workflow
<img width="873" alt="scAEGAN" src="https://user-images.githubusercontent.com/70262340/150944062-c9c72e62-ee8b-41f2-8d97-8d7e8711529a.PNG">

# Summary
scAEGAN is a python based deep learning model that is designed for single-cell-omics and multi-omics integration. scAEGAN does this by using an Autoencoder which learns a low-dimensional embedding of each experiment independently, respecting each sample's uniqueness, protocol. Next, cycleGAN learns a non-linear mapping between these two Autoencoder representations, leveraging the observation that distributions in different latent spaces are similar.


# Installation Requisites 

The required libraries are included in environment.yml file. In order to install these libraries, follow the following steps:

* Creating the conda environment with the following command. This will create and install the libraries included in the environment.yml file for training the scAEGAN.
```
conda env create --prefix ./env --file environment.yml --force
 ```

* The second step is to activate the conda envirnoment. 
```
conda activate ./env      
```


* To evaluate with the given scripts in the Analysis folder on the scAEGAN generated data following libraries are required: 
```
Seurat==4.1.0
clusteval==0.2.1
scclusteval==0.0.0.9
pheatmap==1.0.12
ggplot2==2.3.3.5
cowplot==1.1.1
pdfCluster==1.0-3
cluster==2.1.2
Rtsne==0.15

```
# Usage
There are two steps for the basic usage after activating the conda environment.
*  Training the autoencoder with the given parameters to get the latent representation by running. 
```
python AE.py --input_file <h3 text style="color:#ff0000">Danger</h3> --output_file  --batch_size  --epochs --dropout --learning_rate
```
# Optional arguments
```
--input_file  Specifies the input file to the autoencoder. The input to the autoencoder is in the format of cell by gene matrix. With cells as rows and genes as columns in csv format.
--output_file  Specifies the low dimensional representation of the input from the autoencoder.
--batch_size <int> Specifies the batch size to train the autoencoder. The default batch size is set to 16.
--epochs Specifies <int> the number of epochs for which autoencoder is trained. The default value is set to 200. 
--dropout Specifies the dropout rate used to train the autoencoder. The default value is set to 0.2
--learning_rate Specifies the larning rate. The default learning rate is set to 0.0001
```
*  Training the cyclegan with the given parameters on latent representations obtained from the Autoencoder by running.
   *  The input to the cyclegan is aslo in the format of cell by gene matrix. With cells as rows and genes as columns in csv format

```
python cGANtrain.py --data_path --train_file --test_file --save_path --input_shape --batch_size --epochs
```
# Optional arguments
```
--data_path Specifies the folder path to the training and testing data
--train_file Specifies the training files for training the cGAN for both domains (A and B) that are to be integrated. For instance --train_file train_A.csv train_B.csv
--test_file Specifies the test files. For instance --test_file test_A.csv test_B.csv
--save_path Specifies the folder path where the output from the cGAN in the csv format will be saved.
--input_shape Specifies the shape of the input. The default is set to 50.
--batch_size Specifies the batch size. The default batch size for cGAN is set to 4.
--epochs Specifies the number of epochs for training cGAN. Thed default is set to 200.
```

# Notebooks
* For evaluation, jupyter notebook is given as as example how to run the evaluation on the output from cyclegan.
