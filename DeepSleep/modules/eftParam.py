#                      #  
##                    ##
########################                               
### TTZ/H, Z/H to bb ###
### helper class     ###                               
### for EFT params   ###                               
########################                               
### written by:      ###                               
### Bryan Caraway    ###                               
########################                               
##                    ##                                 
#                      #

import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import os
import re
#import operator as OP
#
import config.ana_cff as cfg
#from makeDatacard import Systematic
#import config.process_norms as p_norms
#from lib.fun_library import weighted_quantile, getZhbbBaseCuts, getZhbbWeight, t2Run
#from lib.TH1 import export1d
#from lib.datacard_shapes import DataCardShapes
#
#import uproot
#import uproot_methods
#from ROOT import TFile, TDirectory, TH1F
#import coffea.hist as hist
from collections import OrderedDict
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt


class EFTFitParams():
    '''
    initialize, calc, and store eft fit 
    of tth/z 0-3 in this class to be accessible
    in EFTParam
    '''
    file_dir = './files/'   # format files/year/mc_files/TT[Z,H]_EFT_val.pkl
    years = cfg.Years
    mc_dir   = 'mc_files/'
    aux_dir = './data/EFT/'
    #
    sig     = ['ttZ','ttH']
    pt_bins = [0,200,300,450,500] # clip at 500
    
    def __init__(self, aux='aux_eftv1.pkl'):
        self.aux_df = pd.read_pickle(f'{self.aux_dir}{aux}')
        self.__worker()

    def __worker(self):
        # assemble eft_df dictionary here by year, signal process, and genZHpt bin
        self.eft_df = {y:{s:{} for s in self.sig} for y in self.years} 
        for y in self.years:
            for s in self.sig:
                if not os.path.exists(f'{self.file_dir}{y}/{self.mc_dir}{s.upper()}_EFT_val.pkl'): continue
                df = pd.read_pickle(f'{self.file_dir}{y}/{self.mc_dir}{s.upper()}_EFT_val.pkl').filter( regex="EFT|genZHpt|NN", axis='columns')
                df['pt_bin'] = pd.cut(df['genZHpt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]-1), bins=self.pt_bins,
                                      labels=[i_bin for i_bin in range(len(self.pt_bins)-1)])
                df = self.calcBeta(df[df['NN']>=0.0])
                # store SM normalized PQR parameters per
                self.eft_df[y][s] = {f'{s}{i}': df[df['pt_bin'] == i].filter(regex=r'c|SM').sum(axis='index')/df[df['pt_bin'] == i]['SM'].sum(axis='index')
                                     for i in range(len(self.pt_bins)-1)}
            #
        #

    def calcBeta(self,df):
         # taken from Jon's code  
         # Build the experiment matrix
         _x =[np.ones(len(self.aux_df.index.values))]
         beta_cols = ['SM']
         for i in range(len(self.aux_df.columns.values)):
             for j in range(i+1):
                 _x.append(self.aux_df.iloc[:,i].values * self.aux_df.iloc[:,j].values)
                 beta_cols.append(f'{self.aux_df.columns.values[i]}_{self.aux_df.columns.values[j]}')
             _x.append(self.aux_df.iloc[:,i].values)
             beta_cols.append(f'{self.aux_df.columns.values[i]}')
         _x = np.matrix(_x).T
         # Build the result matrix y
         _y = np.asmatrix(df.filter(regex=r'EFT').to_numpy()).T
         # Compute beta matrix
         beta = ((_x.T * _x).I * _x.T * _y).A
         return pd.concat([df.reset_index(),pd.DataFrame(data = beta.T, columns=beta_cols)], axis='columns')


class EFTParam():
    '''
    Inheiriting DCParam class to handle the task
    of adding relevent EFT parameter rates
    including WC args, P Q R rates and 
    overall EFT scale factor per gen bin
    
    Must specify ttH or ttZ
    '''
    aux_dir = 'data/EFT'

    def __init__(self):
        #
        self.eft_fit = EFTFitParams().eft_df

    def get_EFT_lines(self, year=None):
        self.wc = None
        dc_lines = []
        for k in self.eft_fit[year]:
            for p,df in self.eft_fit[year][k].items():
                add_lines=[f'# EFT params for process: {p}']
                # create extArg lines for P, Q parameters
                for c in df.keys():
                    if 'SM' in c: continue 
                    add_lines.append(f'{p}_{year}_{c} extArg {np.float32(df[c])}')
                # 
                wc = sorted(set(re.findall(r'c[a-zA-Z]+', ' '.join(df.keys()))))
                self.eft_form(p, year, add_lines, wc)
                # if first iteration, first take care of WC flatParams
                if self.wc is None : 
                    self.wc = wc
                    add_lines = [f'{w} param 0 10' for w in wc] + [f'{w} flatParam' for w in wc] + add_lines
                # add \n to the end of each line
                add_lines = [l+' \n' for l in add_lines]
                dc_lines += add_lines
            #
        #
        return dc_lines
        #

    def eft_form(self, p, y, lines, wc): 
        wc_index = OrderedDict({w: i for i,w in enumerate(wc)})
        o_str = f'CMS_EFT_{p}_{y} rateParam * {p} '
        f_str = ''
        f_index = len(wc_index)
        a_str = ','.join(wc_index.keys())
        for l in lines:
            if '#' in l: continue
            wc2add = re.findall(r'c[a-zA-Z]+',l)
            f_str  += f'@{f_index}'
            a_str  += f",{l.split()[0]}"
            for w in wc2add:
                f_str += f'*@{wc_index[w]}'
            f_str +='+'
            f_index += 1
                        
        #
        o_str += f'({f_str}1.0) {a_str}'
        lines.append(o_str)
        

if __name__ == '__main__':
    eft = EFTParam()
    dc_lines = eft.get_EFT_lines(year='2018')
    print(dc_lines)
    dc_lines = eft.get_EFT_lines(year='2017')
    print('here')
    print(dc_lines)
    
