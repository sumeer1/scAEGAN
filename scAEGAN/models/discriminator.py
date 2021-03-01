from keras.layers import Input
from keras.layers.core import Dense
from keras.models import Model

from .networks_utils import dense_layer


def Discriminator_cGAN(input_shape=(50,), use_wgan=False, use_batch_norm=True, use_leaky_relu=False):
    inputs = Input(shape=input_shape)
    x = inputs

    x = dense_layer(x, units=30, use_batch_norm=use_batch_norm, use_leaky_relu=use_leaky_relu)
    x = dense_layer(x, units=50, use_batch_norm=use_batch_norm, use_leaky_relu=use_leaky_relu)
   # x = dense_layer(x, units=100, use_batch_norm=use_batch_norm, use_leaky_relu=use_leaky_relu)

    activation = None if use_wgan else "sigmoid"
    x = Dense(units=1, activation=activation)(x)

    outputs = x
    return Model(inputs=inputs, outputs=outputs)


def Discriminator(network_type='cGAN', **args):
    assert network_type in {'cGAN'}, "Network'network_type'!!!"

    generators = {
        "cGAN": Discriminator_cGAN,
    }

    return generators[network_type](**args)
