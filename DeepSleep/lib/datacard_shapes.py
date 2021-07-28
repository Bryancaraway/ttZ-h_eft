'''
Define 2d shape histograms 
for datacard per year 
and per datacard channel
'''
import os
import sys
if __name__ == '__main__':
    import subprocess as sb

    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')

import numpy as np
import pandas as pd
import json
import functools
import config.ana_cff as cfg
from lib.fun_library import weighted_quantile, getZhbbWeight, getZhbbBaseCuts, t2Run

class DataCardShapes():
    years = cfg.Years
    file_dir = './files/'
    #
    ref_samples = ['ttZ','ttH']
    #
    nn = cfg.nn
    hist_dict = {}
    #
    def __init__(self, recopt_bins, recosdM_bins, n_NN_bins=10, nn=None, isblind=True): # will have to change isblind for final fit
        self.pt_bins = recopt_bins+[500]
        self.sdM_bins = recosdM_bins
        #self.n_NN_bins = n_NN_bins
        self.nn = cfg.nn if nn is None else nn
        self.isblind = isblind
        #
        self.init_hist_funcs()
    
    @t2Run
    def init_hist_funcs(self):
        for y in  self.years:
            get_pickle= (lambda s: pd.read_pickle(f'{sys.path[1]}/{self.file_dir}{y}/mc_files/{s}_val.pkl'))
            df = pd.concat([get_pickle(s) for s in self.ref_samples], axis='rows', ignore_index=True)
            #df = df[((df[self.nn]>=0.0) & ((df['Hbb']==True) | (df['Zbb'] == True)))]
            df = df[((df[self.nn]>=0.0) & (df['matchedGen_ZHbb_bb']==True))]
            #if y != '2018':
            #    df = df[(df[self.nn] >= 0.0) & (df['process'] != 'ttX')]
            #else:
            #    df2 = pd.read_pickle(f'{self.file_dir}{y}/mc_files/TTZ_bb_val.pkl')
            #    df = pd.concat([df[(df[self.nn] >= 0.0) & (df['process'] == 'ttHbb')],df2[df2[self.nn] >= 0.0]])

            df['genZHpt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]-1, inplace=True)
            df['pt_bin'] = pd.cut(df['genZHpt'], bins=self.pt_bins, 
                                  labels=[f'Zhpt{i_bin}' for i_bin in range(len(self.pt_bins[:-1]))])
            # but we dont really care about pt_bin 0-200, so lets start at index 1
            self.hist_dict[y] = {}
            for i in range(1,len(self.pt_bins[:-1])):
                sub_df = df[df['pt_bin'] == f'Zhpt{i}']
                #quantiles = np.linspace(0,1,self.n_NN_bins+1) # actual quantiles, 10% intervals
                quantiles = [0.0, .05, .25, .35, .50, .70, 1.0] # actual quantiles, 10% intervals

                #if self.isblind:
                nn_df = sub_df[self.nn] 
                #else:
                #    nn_df = sub_df[self.nn][sub_df[self.nn]<=1.7]

                nn_bins = weighted_quantile(nn_df,
                                            quantiles, 
                                            getZhbbWeight(sub_df,y))
                nn_bins[0], nn_bins[-1] = 0,1 # explicity set bin edges to 0 and 1
                if self.isblind == False:
                    #nn_bins = nn_bins[:-1] # up to last bin but mass-side bands
                    nn_bins = nn_bins # dont do anything
                
                #nn_bins = nn_bins[1:] # drop first background dominated bin
                print(nn_bins)
                
                self.hist_dict[y][i] = functools.partial(
                    np.histogram2d,
                    bins=[nn_bins,(self.sdM_bins[i] if type(self.sdM_bins) is dict else self.sdM_bins)])
                #
            #
        #
                                            
    def __getitem__(self,y): # building this like a dictionary 
        try:
            return self.hist_dict[y]
        except KeyError:
            raise KeyError(f"{y} is not a valid top-level key!!!")


class DataCardShapes_NN():
    year = '2018' # just use 2018
    file_dir = './files/'
    ref_samples = cfg.Data_samples
    data_dir = cfg.dnn_ZH_dir

    def __init__(self, n_hl2_bins=64):
        self.out_nn_bins_file = self.data_dir+'/hidden_layer_2_bins.json'
        if os.path.exists(self.out_nn_bins_file):
            self.out_nn_bins_dict = self.get_json_file()
        get_pickle= (lambda s: pd.read_pickle(f'{sys.path[1]}/{self.file_dir}{self.year}/data_files/{s}_val.pkl'))
        self.df = pd.concat([get_pickle(s) for s in self.ref_samples], axis='rows', ignore_index=True)
        self.df = self.df[getZhbbBaseCuts(self.df)].filter(like='NN_')
        #
        self.out_nn_bins_dict = {}
        quantiles = np.linspace(0,1,21)
        for i in range(n_hl2_bins):
            self.out_nn_bins_dict[f'NN_{i}'] = np.quantile(self.df[f'NN_{i}'], quantiles)
            idx = (np.abs(self.out_nn_bins_dict[f'NN_{i}'] - 0.0)).argmin()
            self.out_nn_bins_dict[f'NN_{i}'][idx] = 0.0
            self.out_nn_bins_dict[f'NN_{i}'] = list(np.round(self.out_nn_bins_dict[f'NN_{i}'],3))
        json.dump(self.out_nn_bins_dict, open(self.out_nn_bins_file,'w'), indent=4)

        
    def get_json_file(self):
        _out_dict =  json.load(open(self.out_nn_bins_file,'r'))
        return _out_dict

if __name__ == '__main__':
    test = DataCardShapes([0,200,300,450],[50,80,105,145,200])
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc("figure", figsize=(8, 6*(6./8.)), dpi=200)
    
    print(test['2017'])
    print(test['2017'][2])
    print(test['2017'][2].keywords['bins'])
    #
    fig, ax = plt.subplots()
    fig.subplots_adjust(
        top=0.88,
        bottom=0.11,
        left=0.15,
        right=0.95,
        hspace=0.0,
        wspace=0.0
    )
    y_i = 0.5
    for y in ['2016','2017','2018']:
        for i in [1,2,3]:
            print(f'\n{y} Zhpt bin: {i}')
            x = test[y][i].keywords['bins'][0]
            ax.scatter(x=x,y=y_i*np.ones_like(x),#+test[y][i].keywords['bins'][0],
                        marker='*',
                        label=f'{y}_Zhpt{i}')
            y_i += 0.5
            #
        #
    #
    #ax.legend(bbox_to_anchor=(1.00,1), loc='upper left')
    ax.set_xlim(.1,1.1)
    ax.set_ylim(0,y_i+.5)
    ax.grid(True)
    ax.set_xlabel('DNN bin edge')
    ax.set_ylabel('Datacard Channel')
    ax.set_xticks(np.cumsum([.1]*10))
    ax.set_yticks(np.cumsum([.5]*9))
    ch_dict = {
        1: r'$200 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 300$',
        2: r'$300 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 450$',
        3: r'${p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} > 450$'
    }
    ax.set_yticklabels([f'{y}\n{ch_dict[i]}' for y in ['2016','2017','2018'] for i in [1,2,3]], fontsize='xx-small')
    fig.suptitle('Analysis quantile-defined DNN bins')
    plt.show()
    #sig_df = pd.read_pickle(f'./files/2017/mc_files/TTZH_val.pkl')
    #sub_df = sig_df[(((sig_df['Zh_pt'] >= 300) & (sig_df['Zh_pt']<450)) & ((sig_df['genZHpt'] >= 300)&(sig_df['genZHpt'] < 450)))]
    #hist= test['2017'][2](sub_df[self.nn].to_numpy(), sub_df['Zh_M'].to_numpy(), weights=getZhbbWeight(sub_df,'2017'))
    #print(hist)
