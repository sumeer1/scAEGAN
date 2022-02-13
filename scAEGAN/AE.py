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
import warnings
warnings.filterwarnings('ignore')
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input_file', type=str, required=True)
parser.add_argument('--output_file', type=str, required=True)
parser.add_argument('--dropout_rate', type=int, required=True, default=0.2)
parser.add_argument('--learning_rate', type=int, required=True, default=0.0001)
parser.add_argument('--batch_size', type=int, required=True, default=16)
parser.add_argument('--epochs', type=int, required=True, default=200)
parser.add_argument('--validation_split', type=int, required=False, default=0.2)

args = parser.parse_args()

X = pd.read_csv(args.input_file,sep=',', index_col=0)

model = Sequential()
model.add(Dropout(args.dropout_rate,  input_shape=(X.shape[1],)))
model.add(Dense(300,     activation = 'relu'))
model.add(Dense(50,      activation = 'linear', name = "bottleneck"))
model.add(Dense(300,     activation = 'relu'))
model.add(Dense(X.shape[1],   activation = 'sigmoid'))
model.compile(loss = 'mean_squared_error', optimizer = Adam(lr=args.learning_rate))
model.summary()
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=20)
history = model.fit(X, X, batch_size = args.batch_size, epochs = args.epochs, shuffle = True, verbose = 1, validation_split = args.validation_split, callbacks=[es])
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
RNA_Latent.to_csv(args.output_file,sep=',')
