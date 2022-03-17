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

domainA = pd.read_csv('../DomainA.csv',sep=',', index_col=0)
domainB = pd.read_csv('../DomainB.csv',sep=',', index_col=0)


N_domainA = domainA.values[:,0:(domainA.shape[1])]

N_domainB = domainB.values[:,0:(domainB.shape[1])]

# Input Layer
ncol_domainA = N_domainA.shape[1]
input_dim_domainA = Input(shape = (ncol_domainA, ), name = "domainA")
ncol_domainB = N_domainB.shape[1]
input_dim_domainB = Input(shape = (ncol_domainB, ), name = "domainB")

encoding_dim_domainA = 30
encoding_dim_domainB = 30

dropout_domainA = Dropout(0.2, name = "Dropout_domainA")(input_dim_domainA)
dropout_domainB = Dropout(0.2, name = "Dropout_domainB")(input_dim_scRNAseq2)

encoded_domainA = Dense(encoding_dim_domainA, activation = 'relu', name = "Encoder_domainA")(dropout_domainA)
encoded_domainB = Dense(encoding_dim_domainB, activation = 'relu', name = "Encoder_domainB")(dropout_domainB)

merge = concatenate([encoded_domainA,  encoded_domainB])

bottleneck = Dense(50, kernel_initializer = 'uniform', activation = 'linear', name = "Bottleneck")(merge)

merge_inverse = Dense(encoding_dim_domainA + encoding_dim_domainB, activation = 'relu', name = "Concatenate_Inverse")(bottleneck)

decoded_domainA = Dense(ncol_domainA, activation = 'sigmoid', name = "Decoder_domainA")(merge_inverse)

decoded_domainB = Dense(ncol_domainB, activation = 'sigmoid', name = "Decoder_domainB")(merge_inverse)

autoencoder = Model(input = [input_dim_domainA, input_dim_domainB], output = [decoded_domainA, decoded_domainB])

opt = keras.optimizers.Adam(lr=0.0005)
autoencoder.compile(optimizer = opt, loss={'Decoder_domainA': 'mean_squared_error', 'Decoder_scRNAseq2': 'mean_squared_error', })
autoencoder.summary()

estimator = autoencoder.fit([N_domainA, N_domainB], [N_domainA, N_domainB], epochs = 200, batch_size = 16, validation_split = 0.2, shuffle = True, verbose = 1)
print("Training Loss: ",estimator.history['loss'][-1])
print("Validation Loss: ",estimator.history['val_loss'][-1])
plt.plot(estimator.history['loss'])
plt.plot(estimator.history['val_loss'])
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train','Validation'], loc = 'upper right')
plt.show()

encoder = Model(input = [input_dim_domainA, input_dim_domainB], output = bottleneck)
bottleneck_representation = encoder.predict([N_domainA, N_domainB])

RNA_ATAC_Latent =pd.DataFrame(bottleneck_representation)
RNA_ATAC_Latent.to_csv("../Bottleneck_Representation.csv",sep=',')

