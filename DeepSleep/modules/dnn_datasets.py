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
from lib.fun_library import getZhbbBaseCuts as dnn_cut
from keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder

genmatchreq = 'matchedGen_ZHbb_bb'
#genmatchreq = 'matchedGen_ZHbb_nn'

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
    #dnn_vars = cfg.dnn_ZH_vars 
    #dnn_vars = cfg.nodak8md_dnn_ZH_vars
    #dnn_vars = cfg.withbbvl_dnn_ZH_vars
    dnn_vars = cfg.withbbvl_dnn_ZHgenm_vars
    cut_vars =  ['process','Zh_pt','MET_pt','Zh_M', 'isEleE', 'isMuonE','Zh_bbvLscore', 'passNotHadLep']
    test_train_dir = cfg.dnn_ZH_dir

    def __init__(self):
        self.sig_files = np.array([glob(self.sb_loc+f'{s}_val.pkl') for s in self.sig]).flatten()
        self.bkg_files = np.array([glob(self.sb_loc+f'{b}_val.pkl') for b in self.bkg]).flatten()
        self.s_df , self.b_df = self.get_sigbkg()
        self.prep_class()
        self.sep_test_train()

    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]
        s_df = pd.concat([
            #pd.read_pickle(b).loc[:,pre_vars+['matchedGen_ZHbb']] for b in self.sig_files], ignore_index=True)
            pd.read_pickle(b).loc[:,pre_vars+[genmatchreq]] for b in self.sig_files], ignore_index=True)
        b_df = pd.concat([
            pd.read_pickle(b).loc[:,pre_vars+['tt_type']] for b in self.bkg_files], ignore_index=True)
        #
        #s_df = s_df[s_df['matchedGen_ZHbb'] == True]
        s_df = s_df[s_df[genmatchreq] == True]
        b_df = b_df[((b_df['tt_type'] == 'Semi') & (
            #(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_bb') | (b_df['process'] == 'tt_2b')
            (b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')
        ))]
        return s_df, b_df

    def prep_class(self):
        self.s_df['label'] = 2
        self.b_df['label'] = np.where(self.b_df['process']=='TTBar', 0, 1)
        del self.s_df[genmatchreq], self.b_df['tt_type']
        sb_df = pd.concat([self.s_df,self.b_df])
        encoder = LabelEncoder()
        encoder.fit(sb_df['label'])
        encoded_labels = encoder.transform(sb_df['label'])
        onehot_labels = np_utils.to_categorical(encoded_labels)
        sb_df['label'] = onehot_labels.tolist()
        self.sb_df = sb_df[dnn_cut(sb_df)].sample(frac=1).reset_index(drop=True) # to shuffle dataset
        print(np.unique(self.sb_df['process'],return_counts=True))
        for v in self.cut_vars:
            if v not in self.dnn_vars:
                del self.sb_df[v]

    def sep_test_train(self):
        train_df = self.sb_df.sample(frac=.75, random_state=1)
        test_df  = self.sb_df.drop(train_df.index).copy().sample(frac=1).reset_index(drop=True) 
        train_df.to_pickle(self.test_train_dir+'/trainXY.pkl')
        test_df.to_pickle(self.test_train_dir+'/testXY.pkl')

if __name__ == '__main__':
    _ = DNN_datasets()
