## temp file to store dnn model
import pandas as pd
import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import backend as K

import config.ana_cff as cfg

class DNN_model:

    #alpha  = cfg.dnn_ZH_alpha
    #outDir = cfg.DNNoutputDir
    #outName= cfg.DNNoutputName
    #useW   = cfg.DNNuseWeights
    #

    seq_dict = {
        'Dense'  : (lambda x,y,z: keras.layers.Dense(x, activation='relu', kernel_regularizer=y(z))),
        'Dropout': (lambda x: keras.layers.Dropout(x)),
        'Batch'  : (lambda : keras.layers.BatchNormalization())
    }

    #@classmethod
    #def focal_loss(cls,y_true, y_pred):
    #    gamma = cls.fl_g
    #    alpha = cls.fl_a
    #    pt_1 = tf.where(tf.equal(y_true, 1), y_pred, tf.ones_like(y_pred))
    #    pt_0 = tf.where(tf.equal(y_true, 0), y_pred, tf.zeros_like(y_pred))
    #    return -K.sum(alpha * K.pow(1. - pt_1, gamma) * K.log(pt_1))-K.sum((1-alpha) * K.pow( pt_0, gamma) * K.log(1. - pt_0))
        
    def __init__(self, sequence, other_settings):
        self.sequence = sequence # tuple of ('Dense', (128,keras.regularizers.l1,0.0))
        self.fl_alpha = other_settings['fl_a']
        self.fl_gamma = other_settings['fl_g']
        self.lr_alpha = other_settings['lr_alpha']


    @staticmethod
    def cfl(alpha, gamma=2.):

    alpha = np.array(alpha, dtype=np.float32)

    def categorical_focal_loss_fixed(y_true, y_pred):
        """
        :param y_true: A tensor of the same shape as `y_pred`
        :param y_pred: A tensor resulting from a softmax
        :return: Output tensor.
        """

        # Clip the prediction value to prevent NaN's and Inf's
        epsilon = K.epsilon()
        y_pred = K.clip(y_pred, epsilon, 1. - epsilon)

        # Calculate Cross Entropy
        cross_entropy = -y_true * K.log(y_pred)

        # Calculate Focal Loss
        loss = alpha * K.pow(1 - y_pred, gamma) * cross_entropy

        # Compute mean loss in mini_batch
        return K.mean(K.sum(loss, axis=-1))

    return categorical_focal_loss_fixed


    @staticmethod
    def resetIndex(df_):
        return df_.reset_index(drop=True).copy()


    def Build_Model(self, input_shape ):#mean_,std_):
        #
        main_input = keras.layers.Input(shape=[input_shape], name='input')
        #
        #layer      = keras.layers.Lambda(lambda x: (x - K.constant(mean_)) / K.constant(std_), name='normalizeData')(main_input)
        layer = keras.layers.BatchNormalization()(main_input)
        #
        #layer      = keras.layers.Dense(128, activation='relu', kernel_regularizer=keras.regularizers.l1(0.00))(layer)
        #layer      = keras.layers.Dense(64, activation='relu', kernel_regularizer=keras.regularizers.l1(0.00))(layer)
        #layer      = keras.layers.Dropout(0.5)(layer) 
        for seq in self.sequence:
            layer = self.seq_dict[seq[0]](*seq[1])(layer)
        #
        output     = keras.layers.Dense(3, activation='softmax', name='output')(layer)
            #
        model      = keras.models.Model(inputs=main_input, outputs=output, name='model')
        optimizer  = keras.optimizers.Adam(learning_rate=self.lr_alpha)
        # decay=cfg.dnn_ZH_alpha*1e-3, momentum=0.9, nesterov=True)
        #
        #if (cls.useW and os.path.exists(cls.outDir+cls.outName)): 
        #    model.load_weights(cls.outDir+cls.outName)
        #
        #from focal_loss import BinaryFocalLoss
        model.compile(loss=[self.cfl(self.fl_alpha,self.fl_gamma)], optimizer=optimizer, 
                      metrics=['accuracy',tf.keras.metrics.AUC()])
        ##
        return model

def train_model(m_info):
    m_class = DNN_model(m_info['sequence'],m_info['other_settings'])
    model = m_class.Build_Model(len(cfg.dnn_ZH_vars))
    # --
    # earlystopping here?
    # --
    # open train , val, here?
    # --
    history = model.fit(
        trainX,
        trainY,
        epochs     = m_info['n_epochs'],
        batch_size = m_info['batch_size'],
        validation_data = (valX, valY)
        verbose = 1
    )
    loss ,acc, auc             = model.evaluate(testX, testY)#, sample_weight = testW['DNNweight'].values)
    tr_loss ,tr_acc, tr_auc    = model.evaluate(trainX,trainY)
    val_loss ,val_acc, val_auc = model.evaluate(valX,  valY)
    print("\n\n")
    print(f"Train Set Loss: {tr_loss:10.4f}  Train Set Acc: {tr_acc:10.4f} Train Set AUC: {tr_auc:10.4f}"
    print(f"Val   Set Loss: {val_loss:10.4f}  Val   Set Acc: {val_acc:10.4f} Val   Set AUC: {val_auc:10.4f}"
    print(f"Test  Set Loss: {loss:10.4f}  Test  Set Acc: {acc:10.4f} Test  Set AUC: {auc:10.4f}\n\n"
    


if __name__ == "__main__":
    import json
    import sys
    json_dir = f'sys.path[1]/log/nn/'
    m_info = json.load(open(json_dir/sys.argv[1]))
    train_model(m_info)
