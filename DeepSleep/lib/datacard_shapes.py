'''
Define 2d shape histograms 
for datacard per year 
and per datacard channel
'''
if __name__ == '__main__':
    import subprocess as sb
    import sys
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')

import numpy as np
import pandas as pd
import functools
import config.ana_cff as cfg
from lib.fun_library import weighted_quantile, getZhbbWeight, t2Run

class DataCardShapes():
    years = cfg.Years
    file_dir = './files/'
    #
    ref_sample = 'TTZH_val.pkl'
    #
    hist_dict = {}
    #
    def __init__(self, recopt_bins, recosdM_bins, n_NN_bins=10, isblind=True): # will have to change isblind for final fit
        self.pt_bins = recopt_bins+[500]
        self.sdM_bins = recosdM_bins
        self.n_NN_bins = n_NN_bins
        self.isblind = isblind
        #
        self.init_hist_funcs()
    
    @t2Run
    def init_hist_funcs(self):
        for y in self.years:
            df = pd.read_pickle(f'{self.file_dir}{y}/mc_files/{self.ref_sample}')
            if y != '2018':
                df = df[(df['NN'] >= 0.0) & (df['process'] != 'ttX')]
            else:
                df2 = pd.read_pickle(f'{self.file_dir}{y}/mc_files/TTZ_bb_val.pkl')
                df = pd.concat([df[(df['NN'] >= 0.0) & (df['process'] == 'ttHbb')],df2[df2['NN'] >= 0.0]])

            df['genZHpt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]-1, inplace=True)
            df['pt_bin'] = pd.cut(df['genZHpt'], bins=self.pt_bins, 
                                  labels=[f'Zhpt{i_bin}' for i_bin in range(len(self.pt_bins[:-1]))])
            # but we dont really care about pt_bin 0-200, so lets start at index 1
            self.hist_dict[y] = {}
            for i in range(1,len(self.pt_bins[:-1])):
                sub_df = df[df['pt_bin'] == f'Zhpt{i}']
                quantiles = np.linspace(0,1,self.n_NN_bins+1) # actual quantiles, 10% intervals
                #if self.isblind:
                nn_df = sub_df['NN'] 
                #else:
                #    nn_df = sub_df['NN'][sub_df['NN']<=1.7]

                nn_bins = weighted_quantile(nn_df,
                                            quantiles, 
                                            getZhbbWeight(sub_df,y))
                if self.isblind == False:
                    nn_bins = nn_bins[:6]
                nn_bins = nn_bins[1:] # drop first background dominated bin
                print(nn_bins)
                self.hist_dict[y][i] = functools.partial(
                    np.histogram2d,
                    bins=[nn_bins,self.sdM_bins])
                #
            #
        #
                                            
    def __getitem__(self,y): # building this like a dictionary 
        try:
            return self.hist_dict[y]
        except KeyError:
            raise KeyError(f"{y} is not a valid top-level key!!!")

if __name__ == '__main__':
    test = DataCardShapes([0,200,300,450],[50,80,105,145,200])
    import matplotlib.pyplot as plt
    print(test['2017'])
    print(test['2017'][2])
    y_i = 0
    for y in ['2016','2017','2018']:
        for i in [1,2,3]:
            print(f'\n{y} Zhpt bin: {i}')
            plt.step(x=test[y][i].keywords['bins'][0],y=y_i+test[y][i].keywords['bins'][0],label=f'{y}_Zhpt{i}')
            y_i += 0.5

    plt.legend()
    plt.ylim(0,y_i+2)
    plt.grid(True)
    plt.xlabel('NN bin edge')
    plt.show()
    #sig_df = pd.read_pickle(f'./files/2017/mc_files/TTZH_val.pkl')
    #sub_df = sig_df[(((sig_df['Zh_pt'] >= 300) & (sig_df['Zh_pt']<450)) & ((sig_df['genZHpt'] >= 300)&(sig_df['genZHpt'] < 450)))]
    #hist= test['2017'][2](sub_df['NN'].to_numpy(), sub_df['Zh_M'].to_numpy(), weights=getZhbbWeight(sub_df,'2017'))
    #print(hist)
