from keras.layers import Input
from keras.layers.core import Dense
from keras.models import Model

from .networks_utils import (residual_dense_block,
                             dense_layer)


def Generator_cGAN(input_shape=(50,), use_dropout=True, use_batch_norm=True,
                                 use_leaky_relu=False):
    inputs = Input(shape=input_shape)
    embedding = inputs

    embedding = residual_dense_block(embedding, units=50, use_dropout=use_dropout, use_batch_norm=use_batch_norm,
                                     use_leaky_relu=use_leaky_relu)
    #embedding = residual_dense_block(embedding, units=50, use_dropout=use_dropout, use_batch_norm=use_batch_norm,#use_leaky_relu=use_leaky_relu)

    outputs = Dense(units=input_shape[0], activation=None)(embedding)

    return Model(inputs=inputs, outputs=outputs), inputs, outputs


def Generator(network_type='cGAN', **args):
    assert network_type in {'cGAN'}, "Network 'network_type'!!!"

    generators = {
        "cGAN": Generator_cGAN,
    
    }

    return generators[network_type](**args)
