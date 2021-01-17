## temp file to store dnn model
import pandas as pd
import numpy as np
import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
#import tensorflow as tf
#from tensorflow import keras
#from tensorflow.keras import layers
#from tensorflow.keras import backend as K

from glob import glob
import config.ana_cff as cfg
from keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder

class DNN_datasets:
    '''
    prepare data for 
    training in a 
    multiclassifier
    network
    '''
    sig = ['ttH','ttZ']
    bkg = ['TTBar','ttbb']
    sb_loc = cfg.master_file_path+'/*/mc_files/'
    pre_vars = cfg.dnn_ZH_vars + ['process']

    def __init__(self):
        self.sig_files = np.array([glob(self.sb_loc+f'{s}_val.pkl') for s in self.sig]).flatten()
        self.bkg_files = np.array([glob(self.sb_loc+f'{b}_val.pkl') for b in self.bkg]).flatten()
        self.s_df , self.b_df = self.get_sigbkg()
        self.prep_class()

    def get_sigbkg(self):
        s_df = pd.concat([
            pd.read_pickle(b).loc[:,self.pre_vars+['matchedGenZH']] for b in self.sig_files], ignore_index=True)
        b_df = pd.concat([
            pd.read_pickle(b).loc[:,self.pre_vars+['tt_type']] for b in self.bkg_files], ignore_index=True)
        #
        s_df = s_df[s_df['matchedGenZH'] == True]
        b_df = b_df[((b_df['tt_type'] == 'Semi') & (
            (b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_bb') | (b_df['process'] == 'tt_2b')
        ))]
        print(s_df)
        print(b_df)
        return s_df, b_df

    def prep_class(self):
        self.s_df['label'] = 2
        self.b_df['label'] = np.where(self.b_df['process']=='TTBar', 0, 1)
        del self.s_df['matchedGenZH']
        del self.b_df['process'], self.b_df['tt_type']
        sb_df = pd.concat([self.s_df,self.b_df])
        encoder = LabelEncoder()
        encoder.fit(sb_df['label'])
        encoded_labels = encoder.transform(sb_df['label'])
        onehot_labels = np_utils.to_categorical(encoded_labels)
        print(onehot_labels)
        sb_df['label'] = onehot_labels
        print(sb_df['label'])
    


if __name__ == '__main__':
    _ = DNN_datasets()
