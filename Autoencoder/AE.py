import os
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
from sklearn.manifold import TSNE
from keras.optimizers import Adam
from sklearn.decomposition import PCA
from keras.models import Sequential, Model
from keras import regularizers
import warnings
warnings.filterwarnings('ignore')


X = pd.read_csv('../data.csv',sep=',', index_col=0)

model = Sequential()
model.add(Dropout(0.2,  input_shape=(X.shape[1],)))
model.add(Dense(30,     activation = 'relu'))
model.add(Dense(50,      activation = 'linear', name = "bottleneck"))
model.add(Dense(30,     activation = 'relu'))
model.add(Dense(X.shape[1],   activation = 'sigmoid'))
model.compile(loss = 'mean_squared_error', optimizer = Adam(lr=0.0001))
model.summary()

history = model.fit(X, X, batch_size = 16, epochs = 120, shuffle = True, verbose = 1, validation_split = 0.2)
print("\n" + "Training Accuracy: ", history.history['loss'][-1])
print("Validation Accuracy: ", history.history['val_loss'][-1], "\n")
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Validate'], loc='upper right')
plt.show()

encoder = Model(model.input, model.get_layer('bottleneck').output)
bottleneck_representation = encoder.predict(X)

RNA_Latent =pd.DataFrame(bottleneck_representation)
RNA_Latent.to_csv("../bottleneck_representation.csv",sep=',')
