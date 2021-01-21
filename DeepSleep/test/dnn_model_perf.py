## temp file to store dnn model
import pandas as pd
import numpy as np
import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
from sklearn.feature_selection import f_classif, chi2, mutual_info_classif
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
#import tensorflow as tf
#from tensorflow import keras
#from tensorflow.keras import layers
#from tensorflow.keras import backend as K
from lib.fun_library import save_pdf
import config.ana_cff as cfg

@save_pdf('statcorr_dnn.pdf')
def main():
    trainXY = pd.read_pickle(cfg.dnn_ZH_dir+'/trainXY.pkl')
    fig, ax = plt.subplots()
    trainXY['label'] = np.stack(trainXY['label'].values)[:,2]
    df_target = trainXY['label']
    df        = trainXY
    p_corr = df.corr()
    #sns.set(font_scale=0.001)
    res = sns.heatmap(df.corr(), annot=True, 
                      fmt= '1.2f',annot_kws={"size": 2}, xticklabels=True,
                      yticklabels=True,
                      cmap=plt.cm.Reds, cbar=False, square= False, ax=ax)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 2)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 2)
    fig.suptitle("Pearson Correlation Matrix")
    #plt.show()
    #plt.close()
    #exit()
    #
    #
    chi2_score, chi_2_p_value = chi2(abs(df),df_target)
    f_score, f_p_value = f_classif(df,df_target)
    #mut_info_score = mutual_info_classif(df,df_target)
    #
    stat_sum_df = pd.DataFrame(
        [abs(p_corr['label'].values),abs(chi2_score),abs(chi_2_p_value),abs(f_score),abs(f_p_value)],
        columns = df.keys(), index = ['P_corr','Chi2','Chi2_p','Fscore','Fscore_p']).transpose()
    #stat_sum_df = pd.DataFrame(
    #    [abs(p_corr['label'].values),abs(chi2_score),abs(chi_2_p_value),abs(f_score),abs(f_p_value),abs(mut_info_score)],
    #    columns = df.keys(), index = ['P_corr','Chi2','Chi2_p','Fscore','Fscore_p','mut_info']).transpose()
    plot_heatmap(stat_sum_df, 'Statistical Test Metric Summary', col_norm=False)
    #

def plot_heatmap(score_df,title, norm=False, col_norm=False):
    fig, ax = plt.subplots()
    #plt.figure(figsize=(12,8))
    #sns.set(font_scale=0.1)
    
    res = sns.heatmap(((score_df-score_df.mean())/score_df.std() if col_norm else score_df), 
                      annot=True, fmt= '1.2f',annot_kws={"size": 4}, ax=ax, yticklabels=True,
                      cmap=plt.cm.Reds, cbar=False, square= False, norm=(LogNorm() if norm else None))
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 4)
    fig.suptitle(title)
    #plt.show()
    #plt.close()


if __name__ == '__main__':
    main()
