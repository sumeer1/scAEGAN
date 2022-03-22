"""

author: Sumeer Khan
email:  sameer15khan@gmail.com

"Lower dimension projection of single and multiomics data with autoencoders"

"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from sklearn.manifold import TSNE
from keras.layers import Input, Dense, Dropout
from keras.layers.merge import concatenate
from keras.models import Model
from keras.utils import plot_model
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from keras.optimizers import Adam
from sklearn.decomposition import PCA
from keras.models import Sequential, Model
from keras import regularizers
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
import warnings
warnings.filterwarnings('ignore')
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input_file1', type=str, required=True)
parser.add_argument('--input_file2', type=str, required=True)
parser.add_argument('--output_file1', type=str, required=True)
parser.add_argument('--output_file2', type=str, required=True)
parser.add_argument('--dropout_rate', type=float, required=True, default=0.2)
parser.add_argument('--learning_rate', type=float, required=True, default=0.0001)
parser.add_argument('--batch_size', type=int, required=True, default=16)
parser.add_argument('--epochs', type=int, required=True, default=200)
parser.add_argument('--validation_split', type=float, required=False, default=0.2)

args = parser.parse_args()


################ Mode 1 ########################################################################################################################################################################################

# Input from the domainA (cell by gene matrix)

X = pd.read_csv(args.input_file1,sep=',', index_col=0).transpose()

# Model1
model1 = Sequential()
model1.add(Dropout(args.dropout_rate,  input_shape=(X.shape[1],)))
model1.add(Dense(300,     activation = 'relu'))
model1.add(Dense(50,      activation = 'linear', name = "bottleneck1"))
model1.add(Dense(300,     activation = 'relu'))
model1.add(Dense(X.shape[1],   activation = 'relu'))
model1.compile(loss = 'mean_squared_error', optimizer = Adam(lr=args.learning_rate))
model1.summary()
es1 = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=20)
history1 = model1.fit(X, X, batch_size = args.batch_size, epochs = args.epochs, shuffle = True, verbose = 1, validation_split = args.validation_split,callbacks=[es1]) 
print("\n" + "Training Accuracy: ", history1.history['loss'][-1])
print("Validation Accuracy: ", history1.history['val_loss'][-1], "\n")
plt.plot(history1.history['loss'])
plt.plot(history1.history['val_loss'])
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Validate'], loc='upper right')
plt.show()
#Prediction on domainA
encoder = Model(model1.input, model1.get_layer('bottleneck1').output)
bottleneck_representation1 = encoder.predict(X)
# Saving output of domainA
domainA_Latent =pd.DataFrame(bottleneck_representation1)
domainA_Latent.to_csv(args.output_file1,sep=',')


##################### Model 2 #############################################################################################################################################################################################

# Input from the domainB (cell by gene matrix)

Y = pd.read_csv(args.input_file2,sep=',', index_col=0).transpose()

# Model2
model2 = Sequential()
model2.add(Dropout(args.dropout_rate,  input_shape=(Y.shape[1],)))
model2.add(Dense(300,     activation = 'relu'))
model2.add(Dense(50,      activation = 'linear', name = "bottleneck2"))
model2.add(Dense(300,     activation = 'relu'))
model2.add(Dense(Y.shape[1],   activation = 'relu'))
model2.compile(loss = 'mean_squared_error', optimizer = Adam(lr=args.learning_rate))
model2.summary()
es2 = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=20)
history2 = model1.fit(Y, Y, batch_size = args.batch_size, epochs = args.epochs, shuffle = True, verbose = 1, validation_split = args.validation_split,callbacks=[es2]) 
print("\n" + "Training Accuracy: ", history2.history['loss'][-1])
print("Validation Accuracy: ", history2.history['val_loss'][-1], "\n")
plt.plot(history2.history['loss'])
plt.plot(history2.history['val_loss'])
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Validate'], loc='upper right')
plt.show()
#Prediction on domainB
encoder = Model(model2.input, model2.get_layer('bottleneck2').output)
bottleneck_representation2 = encoder.predict(Y)
# Saving output of domainB
domainB_Latent =pd.DataFrame(bottleneck_representation2)
domainB_Latent.to_csv(args.output_file2,sep=',')










































encoder = Model(model1.input, model1.get_layer('bottleneck1').output)
bottleneck_representation1 = encoder.predict(X)

RNA_Latent =pd.DataFrame(bottleneck_representation1)
RNA_Latent.to_csv(args.output_file1,sep=',')

