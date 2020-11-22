import os
import keras
import numpy as np
import pandas as pd
from umap import UMAP
import matplotlib as mpl
from sklearn.manifold import TSNE
from keras.layers import Input, Dense, Dropout
from keras.layers.merge import concatenate
from keras.models import Model
from keras.utils import plot_model
import matplotlib.pyplot as plt
from numpy.random import seed
seed(1)
from tensorflow import set_random_seed
set_random_seed(2)
import warnings
warnings.filterwarnings('ignore')

scRNAseq1 = pd.read_csv('../DomainA.csv',sep=',', index_col=0).transpose()
scRNAseq2 = pd.read_csv('../DomainB.csv',sep=',', index_col=0).transpose()


X_scRNAseq1 = scRNAseq1.values[:,0:(scRNAseq1.shape[1])]

X_scRNAseq2 = scRNAseq2.values[:,0:(scRNAseq2.shape[1])]

# Input Layer
ncol_scRNAseq1 = X_scRNAseq1.shape[1]
input_dim_scRNAseq1 = Input(shape = (ncol_scRNAseq1, ), name = "scRNAseq1")
ncol_scRNAseq2 = scRNAseq2.shape[1]
input_dim_scRNAseq2 = Input(shape = (ncol_scRNAseq2, ), name = "scRNAseq2")

encoding_dim_scRNAseq1 = 30
encoding_dim_scRNAseq2 = 30

dropout_scRNAseq1 = Dropout(0.2, name = "Dropout_scRNAseq1")(input_dim_scRNAseq1)
dropout_scRNAseq2 = Dropout(0.2, name = "Dropout_scRNAseq2")(input_dim_scRNAseq2)

encoded_scRNAseq1 = Dense(encoding_dim_scRNAseq1, activation = 'relu', name = "Encoder_scRNAseq1")(dropout_scRNAseq1)
encoded_scRNAseq2 = Dense(encoding_dim_scRNAseq2, activation = 'relu', name = "Encoder_scRNAseq2")(dropout_scRNAseq2)

merge = concatenate([encoded_scRNAseq1,  encoded_scRNAseq2])

bottleneck = Dense(50, kernel_initializer = 'uniform', activation = 'linear', name = "Bottleneck")(merge)

merge_inverse = Dense(encoding_dim_scRNAseq1 + encoding_dim_scRNAseq2, activation = 'relu', name = "Concatenate_Inverse")(bottleneck)

decoded_scRNAseq1 = Dense(ncol_scRNAseq1, activation = 'sigmoid', name = "Decoder_scRNAseq1")(merge_inverse)

decoded_scRNAseq2 = Dense(ncol_scRNAseq2, activation = 'sigmoid', name = "Decoder_scRNAseq2")(merge_inverse)

autoencoder = Model(input = [input_dim_scRNAseq1, input_dim_scRNAseq2], output = [decoded_scRNAseq1, decoded_scRNAseq2])

opt = keras.optimizers.Adam(lr=0.0005)
autoencoder.compile(optimizer = opt, loss={'Decoder_scRNAseq1': 'mean_squared_error', 'Decoder_scRNAseq2': 'mean_squared_error', })
autoencoder.summary()

plot_model(autoencoder, to_file='autoencoder_graph.png')

estimator = autoencoder.fit([X_scRNAseq1, X_scRNAseq2], [X_scRNAseq1, X_scRNAseq2], epochs = 200, batch_size = 16, validation_split = 0.2, shuffle = True, verbose = 1)
print("Training Loss: ",estimator.history['loss'][-1])
print("Validation Loss: ",estimator.history['val_loss'][-1])
plt.plot(estimator.history['loss'])
plt.plot(estimator.history['val_loss'])
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train','Validation'], loc = 'upper right')
plt.show()

encoder = Model(input = [input_dim_scRNAseq1, input_dim_scRNAseq2], output = bottleneck)
bottleneck_representation = encoder.predict([X_scRNAseq1, scRNAseq2])
print(pd.DataFrame(bottleneck_representation).shape)
print(pd.DataFrame(bottleneck_representation).iloc[0:5,0:5])

RNA_ATAC_Latent =pd.DataFrame(bottleneck_representation)
RNA_ATAC_Latent.to_csv("../Bottleneck_Representation.csv",sep=',')

