## temp file to store dnn model
import pandas as pd
import os

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import backend as K


import cfg.deepsleepcfg as cfg

class DNN_model:

    alpha  = cfg.dnn_ZH_alpha
    outDir = cfg.DNNoutputDir
    outName= cfg.DNNoutputName
    useW   = cfg.DNNuseWeights
    #
    fl_g = cfg.fl_gamma
    fl_a = cfg.fl_alpha

    @classmethod
    def focal_loss(cls,y_true, y_pred):
        gamma = cls.fl_g
        alpha = cls.fl_a
        pt_1 = tf.where(tf.equal(y_true, 1), y_pred, tf.ones_like(y_pred))
        pt_0 = tf.where(tf.equal(y_true, 0), y_pred, tf.zeros_like(y_pred))
        return -K.sum(alpha * K.pow(1. - pt_1, gamma) * K.log(pt_1))-K.sum((1-alpha) * K.pow( pt_0, gamma) * K.log(1. - pt_0))
        
    @staticmethod
    def resetIndex(df_):
        return df_.reset_index(drop=True).copy()

    @classmethod
    def Build_Model(cls, inputs_,outputs_,mean_,std_):
        #
        main_input = keras.layers.Input(shape=[inputs_], name='input')
        #
        layer      = keras.layers.Lambda(lambda x: (x - K.constant(mean_)) / K.constant(std_), name='normalizeData')(main_input)
        #
        layer      = keras.layers.Dense(128, activation='relu', kernel_regularizer=keras.regularizers.l1(0.00))(layer)
        layer      = keras.layers.Dense(64, activation='relu', kernel_regularizer=keras.regularizers.l1(0.00))(layer)
        layer      = keras.layers.Dropout(0.5)(layer) 
        #
        output     = keras.layers.Dense(1, activation='sigmoid', name='output')(layer)
            #
        model      = keras.models.Model(inputs=main_input, outputs=output, name='model')
        optimizer  = keras.optimizers.Adam(learning_rate=cls.alpha)
        # decay=cfg.dnn_ZH_alpha*1e-3, momentum=0.9, nesterov=True)
        #
        if (cls.useW and os.path.exists(cls.outDir+cls.outName)): 
            model.load_weights(cls.outDir+cls.outName)
        #
        #from focal_loss import BinaryFocalLoss
        model.compile(loss=[cls.focal_loss], optimizer=optimizer, 
                      metrics=['accuracy',tf.keras.metrics.AUC()])
        ##
        return model
