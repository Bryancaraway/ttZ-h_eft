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
    #trainXY = pd.read_pickle(cfg.dnn_ZH_dir+'/oldtrainXY.pkl')
    trainXY = pd.read_pickle(cfg.dnn_ZH_dir+'/trainXY.pkl')
    trainXY['sig_label'] = np.stack(trainXY['label'].values)[:,2]
    trainXY['ttbar_label'] = np.stack(trainXY['label'].values)[:,0]
    trainXY['ttbb_label'] = np.stack(trainXY['label'].values)[:,1]
    del trainXY['label']
    print(sum(trainXY['sig_label']), sum(trainXY['ttbb_label']))


    def test_target(target_label, title, _df):
        df_target = _df[target_label]
        df        = _df
        p_corr = df.corr()
        #_out = pd.DataFrame(p_corr[target_label].values, index=df.keys(), columns=[f'{title}_Corr'])
        #_out.drop(index=['sig_label','ttbb_label','ttbar_label'] ,inplace=True)
        #_out[f'{title}_rank'] = abs(_out[f'{title}_Corr']).rank() 
        #return _out
        sns.set(font_scale=0.001)
        fig, ax = plt.subplots()
        res = sns.heatmap(df.corr(), annot=True, 
                          fmt= '1.2f',annot_kws={"size": 2}, xticklabels=True,
                          yticklabels=True,
                          cmap=plt.cm.Reds, cbar=False, square= False, ax=ax)
        res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 2)
        res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 2)
        fig.suptitle(f"Pearson Correlation Matrix ({title})")
        #
        #
        #chi2_score, chi_2_p_value = chi2(abs(df),df_target)
        #_out = pd.DataFrame(chi2_score, index=df.keys(), columns=[f'{title} Mut Info'])
        #return _out

        #f_score, f_p_value = f_classif(df,df_target)        
        mut_info_score = mutual_info_classif(df,df_target)
        _out = pd.DataFrame(mut_info_score, index=df.keys(), columns=[f'{title}_Mut_Info'])
        _out.drop(index=['sig_label','ttbb_label','ttbar_label'] ,inplace=True)
        _out[f'{title}_rank'] = _out[f'{title}_Mut_Info'].rank() 
        return _out
        #stat_sum_df = pd.DataFrame(
        #    [abs(p_corr[target_label].values),abs(chi2_score),abs(chi_2_p_value),abs(f_score),abs(f_p_value)],
        #    columns = df.keys(), index = ['P_corr','Chi2','Chi2_p','Fscore','Fscore_p']).transpose()
        #stat_sum_df = pd.DataFrame(
        #    [abs(p_corr[target_label].values),abs(chi2_score),10**(chi_2_p_value),abs(f_score),10**(f_p_value),abs(mut_info_score)],
        #    columns = df.keys(), index = ['P_corr','Chi2','10^_Chi2_p','Fscore','10^_Fscore_p','mut_info']).transpose()
        #plot_heatmap(stat_sum_df, f'Statistical Test Metric Summary ({title})', col_norm=False)
    #
    # compare sig vs all
    testing_df = pd.concat([trainXY[trainXY['sig_label'] == 0].sample(sum(trainXY['sig_label'] == 1)), trainXY[trainXY['sig_label'] == 1]],ignore_index=True) # equal amounts of signal and background
    mut_info = test_target('sig_label', 'SIGvsBKG', testing_df)
    #test_target('sig_label', 'SIGvsBKG', testing_df)
    #print(mut_info)
    del testing_df
    # now compare only looking at sig vs tt+bb
    testing_df = pd.concat([trainXY[trainXY['ttbb_label'] == 1], trainXY[trainXY['sig_label'] == 1].sample(sum(trainXY['ttbb_label'] == 1))],ignore_index=True) # equal amounts of signal and tt+bb
    mut_info = pd.concat([mut_info,test_target('sig_label', 'SIGvsTTBB', testing_df)], axis='columns')
    #test_target('sig_label', 'SIGvsTTBB', testing_df)
    #print(mut_info)
    del testing_df
    # now compare only looking at sig vs ttbar
    testing_df = pd.concat([trainXY[trainXY['ttbar_label'] == 1], trainXY[trainXY['sig_label'] == 1].sample(sum(trainXY['ttbar_label'] == 1))],ignore_index=True) # equal amounts of signal and ttbar
    mut_info = pd.concat([mut_info, test_target('sig_label', 'SIGvsTTBar', testing_df)], axis='columns')
    #test_target('sig_label', 'SIGvsTTBar', testing_df)
    #print(mut_info)
    del testing_df
    #
    mut_info.sort_values(by='SIGvsTTBB_rank', inplace=True)
    mut_info.index.name = 'Inputs'
    #plot_heatmap(mut_info, 'Pear Corr', cmap=False)
    plot_heatmap(mut_info, 'Mut Info', cmap=False)

            

def plot_heatmap(score_df,title, norm=False, col_norm=False, cmap=True):
    fig, ax = plt.subplots()
    #plt.figure(figsize=(12,8))
    #sns.set(font_scale=0.1)
    from matplotlib.colors import ListedColormap
    res = sns.heatmap(((score_df-score_df.mean())/score_df.std() if col_norm else score_df), 
                      annot=True, fmt= '1.4f',annot_kws={"size": 4}, ax=ax, yticklabels=True,
                      linewidths = None if cmap else .05, linecolor= None if cmap else 'k' ,
                      cmap=plt.cm.Reds if cmap else ListedColormap(['white']), 
                      cbar=False, square= False, norm=(LogNorm() if norm else None))
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 4)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 4, rotation=0)
    fig.suptitle(title)
    #plt.show()
    #plt.close()


if __name__ == '__main__':
    main()
