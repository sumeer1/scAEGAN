"""

author: Sumeer Khan
email:  sameer15khan@gmail.com

"Lower dimension projection of single and multiomics data with autoencoders"

"""

import os
os.environ['KERAS_BACKEND'] = 'tensorflow'
import numpy as np
import pandas as pd
import matplotlib as mpl
from keras.layers import Input, Dense, Dropout
from keras.layers.merge import concatenate
from keras import backend as K
from keras.models import Model
from keras.utils import plot_model
import matplotlib.pyplot as plt
from keras.optimizers import Adam
from keras.constraints import UnitNorm, Constraint
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


################ Model 1 ########################################################################################################################################################################################

# Input from the domainA (cell by gene matrix)

X = pd.read_csv(args.input_file1,sep=',', index_col=0).transpose()


class WeightsOrthogonalityConstraint(Constraint):
    def __init__(self, encoding_dim, weightage = 1.0, axis = 0):
        self.encoding_dim = encoding_dim
        self.weightage = weightage
        self.axis = axis
        
    def weights_orthogonality(self, w):
        if(self.axis==1):
            w = K.transpose(w)
        if(self.encoding_dim > 1):
            m = K.dot(K.transpose(w), w) - K.eye(self.encoding_dim)
            return self.weightage * K.sqrt(K.sum(K.square(m)))
        else:
            m = K.sum(w ** 2) - 1.
            return m

    def __call__(self, w):
        return self.weights_orthogonality(w)


# Model1
model1 = Sequential()
model1.add(Dropout(args.dropout_rate,  input_shape=(X.shape[1],)))
model1.add(Dense(300,     activation = 'relu', use_bias=True, kernel_regularizer=WeightsOrthogonalityConstraint(300, weightage=1., axis=0))) 
model1.add(Dense(50,      activation = 'linear', name = "bottleneck1"))
model1.add(Dense(300,     activation = 'relu'))
model1.add(Dense(X.shape[1],   activation = 'sigmoid'))
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
encoder1 = Model(model1.input, model1.get_layer('bottleneck1').output)
bottleneck_representation1 = encoder1.predict(X)
# Saving output of domainA
domainA_Latent =pd.DataFrame(bottleneck_representation1)
domainA_Latent.to_csv(args.output_file1,sep=',')


##################### Model 2 #############################################################################################################################################################################################

# Input from the domainB (cell by gene matrix)

Y = pd.read_csv(args.input_file2,sep=',', index_col=0).transpose()

# Model2
model2 = Sequential()
model2.add(Dropout(args.dropout_rate,  input_shape=(Y.shape[1],)))
model2.add(Dense(300,     activation = 'relu', use_bias=True, kernel_regularizer=WeightsOrthogonalityConstraint(300, weightage=1., axis=0)))
model2.add(Dense(50,      activation = 'linear', name = "bottleneck2"))
model2.add(Dense(300,     activation = 'relu'))
model2.add(Dense(Y.shape[1],   activation = 'sigmoid'))
model2.compile(loss = 'mean_squared_error', optimizer = Adam(lr=args.learning_rate))
model2.summary()
es2 = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=20)
history2 = model2.fit(Y, Y, batch_size = args.batch_size, epochs = args.epochs, shuffle = True, verbose = 1, validation_split = args.validation_split,callbacks=[es2]) 
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
encoder2 = Model(model2.input, model2.get_layer('bottleneck2').output)
bottleneck_representation2 = encoder2.predict(Y)
# Saving output of domainB
domainB_Latent =pd.DataFrame(bottleneck_representation2)
domainB_Latent.to_csv(args.output_file2,sep=',')

