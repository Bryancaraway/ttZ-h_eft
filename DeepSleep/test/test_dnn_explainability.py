## temp file to store dnn model
import pandas as pd
import numpy as np
import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')

from sklearn.feature_selection import f_classif, chi2, mutual_info_classif
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
from scipy.cluster import hierarchy


import matplotlib.pyplot as plt
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
#import tensorflow as tf
#from tensorflow import keras
#from tensorflow.keras import layers
#from tensorflow.keras import backend as K
from lib.fun_library import save_pdf, getZhbbBaseCuts
import config.ana_cff as cfg
from modules.dnn_model import DNN_model
#
import multiprocessing
from functools import partial

class DNN_Explainability():
    def __init__(self):
        self.getModel() # build model, load weights
        self.prepData() # get signal test data
        #
        save_pdf('dnn_input_hierarchy.pdf')(self.plot_input_hierarchy)()
        #weights = self.model.get_weights()
        #for i in range(len(weights[:-2])):
        #    print('\n')
        #    print(weights[i])
        #    print(weights[i].shape)
        #    print('\n')
        pred = self.model.predict(self.sigX)[:,2]
        #print(self.new_model.predict(self.sigX.iloc[0:5]))
        plt.hist(self.new_model.predict(self.sigX)[:,0])
        plt.show()

        #idx = 4
        #cut = (pred>0.9) & (self.sigY == 0)
        #from lime import lime_tabular
        #explainer = lime_tabular.LimeTabularExplainer(
        #    self.sigX.to_numpy(), mode="classification",
        #    #class_names=['tt+LF,tt+cc','tt+bb','wm Signal'],
        #    feature_names=self.sigX.columns,
        #    verbose=True,
        #)
        #print(self.sigX[cut].iloc[idx])
        #print(pred[cut][idx])
        #print(self.sigY[cut][idx])
        #explanation = explainer.explain_instance(self.sigX[cut].iloc[idx], self.model.predict, top_labels=3, num_features=20,)
        ##print(explanation.available_labels())
        #
        #explanation.as_pyplot_figure(2)
        #plt.tight_layout()
        #plt.show()
        #exit()
        auc = roc_auc_score(self.sigY, pred)
        print(auc)
        self.corrupted_pred = {}
        self.corrupted_auc = {}
        # 
        pool = multiprocessing.Pool()
        #
        for var in self.sigX.columns:
            self.worker(var)
        self.get_max_corrupt()
        #_ = map(self.worker, self.sigX.columns)
        # calculate permutation importance
        
        #print("Permutation importance of well-matched signal and tt+bb, tt+LF, tt+cc events")
        print("Permutation importance, shuffling for background only")
        print("Prediction - Corrupted Prediction: (mean) +/- (std)")
        # sort dict by value
        for var in self.corrupted_pred:
            mean, std = np.mean(pred-self.corrupted_pred[var]), np.std(pred-self.corrupted_pred[var])
            print(f'{var:<25}{mean:.3f} +/- {std:.3f}')
            self.corrupted_auc[var] = roc_auc_score(self.sigY, self.corrupted_pred[var])
        self.corrupted_auc = {k:v for k,v in sorted(self.corrupted_auc.items(), key=(lambda _: abs(auc - _[1])))}
        print("\n")
        print("AUC - Corrupted AUC, and % difference")
        for var in self.corrupted_auc:
            corr_auc = self.corrupted_auc[var]
            print(f'{var:<25}{auc-corr_auc:.3f},    {100*(auc-corr_auc)/auc:.3f} %')
    

    
    def getModel(self):
        dnn_weights = 'withbbvlnewgenm_model.h5'
        dnn_vars   = cfg.withbbvl_dnn_ZHgenm_vars
        m_info =  {
            'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 
            'other_settings': {'fl_a': [0.75, 1, 0.25], 'fl_g': 0.25, 'lr_alpha': 0.0003}, 
            'n_epochs': 100, 'batch_size': 10256,
        }
        dnn = DNN_model(m_info['sequence'],m_info['other_settings']) 
        self.model = dnn.Build_Model(len(dnn_vars), load_weights=dnn_weights)#'nn_ttzh_model.h5')    
        self.new_model = dnn.Build_New_Model(len(dnn_vars), self.model)


    def prepData(self):
        testXY = pd.read_pickle(cfg.dnn_ZH_dir+'/oldtestXY.pkl')
        resetIndex = (lambda df: df.reset_index(drop=True).copy())
        testX ,testY = testXY.drop(columns=['label']), np.stack(resetIndex(testXY['label']).values)[:,2]
        #self.sigX = resetIndex(testX[testY == 1])
        self.sigX = resetIndex(testX)
        #self.sigY = testY[testY == 1]
        self.sigY = testY
        #print(testX.columns)


    def get_max_corrupt(self):
        corrupted_data = self.sigX.copy()
        for i,var in enumerate(self.sigX.columns):
            #corrupted_data.loc[:,var] = corrupted_data.loc[:,var].sample(frac=1,random_state=i).values
            corrupted_data.loc[(self.sigY==1),var] = corrupted_data.loc[(self.sigY==1),var].sample(frac=1,random_state=i).values
        self.corrupted_pred['max'] = np.array(self.model.predict(corrupted_data)[:,2])
            


    def worker(self, var):
        _out = None
        n_perm = 1
        for i in range(n_perm): # run over 10  random shuffles
            corrupted_data = self.sigX.copy()
            #corrupted_data.loc[:,var] = corrupted_data[var].sample(frac=1,random_state=i).values
            # corrupted according to true class
            #corrupted_data.loc[(self.sigY==1),var] = corrupted_data.loc[(self.sigY==1),var].sample(frac=1,random_state=i).values
            corrupted_data.loc[(self.sigY==0),var] = corrupted_data.loc[(self.sigY==0),var].sample(frac=1,random_state=i).values

            #if var in corrupted_pred:
            if _out is not None:
                #corrupted_pred[var] += np.array(model.predict(corrupted_data)[:,2])
                _out += np.array(self.model.predict(corrupted_data)[:,2])
            else:
                #corrupted_pred[var] = np.array(model.predict(corrupted_data)[:,2])
                _out = np.array(self.model.predict(corrupted_data)[:,2])
        self.corrupted_pred[var] = _out/n_perm
        #_out = _out/n_perm
        #return {var:_out}


    def plot_input_hierarchy(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
        fig.subplots_adjust( bottom=0.20)

        corr = spearmanr(self.sigX.to_numpy()).correlation
        corr_linkage = hierarchy.ward(corr)
        dendro = hierarchy.dendrogram(
            corr_linkage, labels=self.sigX.columns, ax=ax1, leaf_rotation=90,
        )
        dendro_idx = np.arange(0, len(dendro['ivl']))
        

        for l in ax1.get_xticklabels():
            l.set_fontsize(6)
        ax2.imshow(corr[dendro['leaves'], :][:, dendro['leaves']])
        ax2.set_xticks(dendro_idx)
        ax2.set_yticks(dendro_idx)
        ax2.set_xticklabels(dendro['ivl'], rotation='vertical',fontsize=6)
        ax2.set_yticklabels(dendro['ivl'],fontsize=6)
        fig.tight_layout()
        from collections import defaultdict
        def print_hier_thresh(thresh=1):
            reduced_idx = hierarchy.fcluster(corr_linkage, thresh, criterion='distance')
            idx_mapping = defaultdict(list)
            for i_, id_ in enumerate(reduced_idx):
                idx_mapping[id_].append(i_)
            indexes = [v[0] for v in idx_mapping.values()]
            print(f"reduced inputs: {len(indexes)}, with thresh: {thresh}\n", self.sigX.columns[indexes])
            return indexes
        #
        vars1p0  = print_hier_thresh(1)
        vars1p25 = print_hier_thresh(1.25)
        vars1p5  = print_hier_thresh(1.5)
        print("extra vars in 1p25 not in 1p5 :\n",self.sigX.columns[[v for v in vars1p25 if v not in vars1p5]])
        print("extra vars in 1p0  not in 1p25:\n",self.sigX.columns[[v for v in vars1p0 if v not in vars1p25]])
        print("extra vars in all  not in 1p0 :\n",self.sigX.columns[[v for v in range(len(self.sigX.columns)) if v not in vars1p0]])

if __name__ == '__main__':
    DNN_Explainability()
