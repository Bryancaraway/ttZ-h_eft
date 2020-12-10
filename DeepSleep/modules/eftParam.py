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
    
    def __init__(self):
        # seperate aux for different signals
        self.aux_df = {s : pd.read_pickle(f'{self.aux_dir}/aux_eft_{s.upper()}.pkl') for s in self.sig }
        self.__worker()

    def __worker(self):
        # assemble eft_df dictionary here by year, signal process, and genZHpt bin
        self.eft_df = {y:{s:{} for s in self.sig} for y in self.years} 
        cut_for_fit = (lambda x: ((x['NN']>=0.0) & (x['EFT183'] < 100)) )
        for y in self.years:
            for s in self.sig:
                if not os.path.exists(f'{self.file_dir}{y}/{self.mc_dir}{s.upper()}_EFT_val.pkl'): continue
                df = pd.read_pickle(f'{self.file_dir}{y}/{self.mc_dir}{s.upper()}_EFT_val.pkl').filter( regex="EFT|genZHpt|NN", axis='columns')
                df['pt_bin'] = pd.cut(df['genZHpt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]-1), bins=self.pt_bins,
                                      labels=[i_bin for i_bin in range(len(self.pt_bins)-1)])
                #df = self.calcBeta(df[df['NN']>=0.0], s)
                df = self.calcBeta(df[cut_for_fit(df)], s)
                # store SM normalized PQR parameters per
                self.eft_df[y][s] = {f'{s}{i}': df[df['pt_bin'] == i].filter(regex=r'c|SM').sum(axis='index')/df[df['pt_bin'] == i]['SM'].sum(axis='index')
                                     for i in range(len(self.pt_bins)-1)}
            #
        #

    def calcBeta(self,df, sample):
        # taken from Jon's code  
        # Build the experiment matrix
        aux_df = self.aux_df[sample]
        _x =[np.ones(len(aux_df.index.values))]
        beta_cols = ['SM']
        for i in range(len(aux_df.columns.values)):
            for j in range(i+1):
                _x.append(aux_df.iloc[:,i].values * aux_df.iloc[:,j].values)
                beta_cols.append(f'{aux_df.columns.values[i]}_{aux_df.columns.values[j]}')
            _x.append(aux_df.iloc[:,i].values)
            beta_cols.append(f'{aux_df.columns.values[i]}')
        _x = np.matrix(_x).T
        # Build the result matrix y
        _y = np.asmatrix(df.filter(regex=r'EFT').to_numpy()).T
        # Compute beta matrix
        beta = ((_x.T * _x).I * _x.T * _y).A
        return pd.concat([df.reset_index(),pd.DataFrame(data = beta.T, columns=beta_cols)], axis='columns')

class BkgEFTFitParams(EFTFitParams):
    # copy the functionality from EFTFItParams
    bkg = ['ttbb']
    aux_dict = {'ttbb':'aux_eft_TTZ.pkl'}

    def __init__(self):
        self.aux_df = {s : pd.read_pickle(f'{self.aux_dir}/{self.aux_dict[s]}') for s in self.bkg}
        self.__worker()

    def __worker(self, kinem=['Zh_pt','ttbb_genbb_pt','ttbb_genbb_invm','process'], k_bins=[-np.inf,np.inf]):
        # assemble eft_df dictionary here by year, signal process, and genZHpt bin
        self.eft_df = {y:{s:{} for s in self.bkg} for y in self.years} 
        cut_for_fit = (lambda x: ((x['NN']>=0.0) & (x['EFT183'] < 100)) )
        for y in self.years:
            for s in self.bkg:
                if not os.path.exists(f'{self.file_dir}{y}/{self.mc_dir}{s.upper()}_EFT_val.pkl'): continue
                df = pd.read_pickle(f'{self.file_dir}{y}/{self.mc_dir}{s.upper()}_EFT_val.pkl').filter( regex=f"EFT|{'|'.join(kinem)}|NN", axis='columns')
                #df['pt_bin'] = pd.cut(df['genZHpt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]-1), bins=self.pt_bins,
                #labels=[i_bin for i_bin in range(len(self.pt_bins)-1)])
                #df = self.calcBeta(df[df['NN']>=0.0], s)
                df = self.calcBeta(df[cut_for_fit(df)], s)
                # store SM normalized PQR parameters per
                #self.eft_df[y][s] = {f'{s}{i}': df[df['pt_bin'] == i].filter(regex=r'c|SM').sum(axis='index')/df[df['pt_bin'] == i]['SM'].sum(axis='index')
                #                     for i in range(len(self.pt_bins)-1)}
                self.eft_df[y][s] = df.filter(regex=r'^(?!EFT)')
            #
        #

        

class EFTParam():
    '''
    Inheiriting DCParam class to handle the task
    of adding relevent EFT parameter rates
    including WC args, P Q R rates and 
    overall EFT scale factor per gen bin
    
    Must specify ttH or ttZ
    '''
    aux_dir = 'data/EFT'
    out_dict = {}
    wc_range = {'cQei' :[-200,200],
                'cQl3i':[-200,200],
                'cQlMi':[-200,200],
                'cbW'  :[-10,10],
                'cpQ3' :[-20,22],
                'cpQM' :[-30,50],
                'cpt'  :[-40,40],
                'cptb' :[-40,40],
                'ctG'  :[-3,3],
                'ctW'  :[-6,6],
                'ctZ'  :[-7,7],
                'ctei' :[-200,200],
                'ctlSi':[-200,200],
                'ctlTi':[-200,200],
                'ctli' :[-200,200],
                'ctp'  :[-35,65]}


    def __init__(self):
        #
        self.eft_fit = EFTFitParams().eft_df
        self.eftbkg_fit  = BkgEFTFitParams().eft_df
        #print(BkgEFTFitParams().eft_df['2018']['ttbb'].keys())

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
                wc = sorted(set(re.findall(r'c[a-zA-Z3]+', ' '.join(df.keys()))))
                self.eft_formula(p, year, add_lines, wc)
                # if first iteration, first take care of WC flatParams
                if self.wc is None : 
                    self.wc = wc
                    add_lines = [f'{w} param 0 {self.wc_range[w][1]}' for w in wc] + [f'{w} flatParam' for w in wc] + add_lines
                # add \n to the end of each line
                add_lines = [l+' \n' for l in add_lines]
                dc_lines += add_lines
            #
        #
        return dc_lines
        #

    def eft_formula(self, p, y, lines, wc): 
        wc_index = OrderedDict({w: i for i,w in enumerate(wc)})
        o_str = f'CMS_EFT_{p}_{y} rateParam * {p} '
        f_str = ''
        f_index = len(wc_index)
        a_str = ','.join(wc_index.keys())
        for l in lines:
            if '#' in l: continue
            f_str  += f'{l.split()[-1]}'
            # special format for parameter values
            #pq_val = l.split()[-1]
            #f_str += re.search('\d+\.\d{5}',pq_val).group() + (re.search(r'e.*',pq_val).group() if 'e' in pq_val else '')
            #
            #f_str  += f'@{f_index}'
            #a_str  += f",{l.split()[-1]}"
            #a_str  += f",{l.split()[0]}"
            wc2add = re.findall(r'c[a-zA-Z3]+',l)
            for w in wc2add:
                f_str += f'*@{wc_index[w]}'
            f_str +='+'
            f_index += 1
                        
        #
        o_str += f'({f_str}1.0) {a_str}'
        lines.append(o_str)

    def save_to_dict(self, year=None, pt_bins=[1,2,3], force_year=None): # construct dictionary based on datacard channels
        if force_year is None: force_year = year
        self.wc = None
        out_dict  = {}
        for k in self.eft_fit[force_year]: # (ttZ,ttH)
            for p,df in self.eft_fit[force_year][k].items(): # (ttZ0,ttZ1,ttZ2,ttZ3)
                wc_dict = {}
                for wc in df.keys():
                    pqr = wc.split('_')
                    if len(pqr) > 1:
                        wc_dict[tuple(pqr)] = df[wc]
                    elif 'SM' in wc:
                        wc_dict[('sm','sm')] = df[wc]
                    else:
                        wc_dict[('sm',wc)] = df[wc]
                for pt in pt_bins:
                    out_dict[(p,f'y{year}_Zhpt{pt}')] = wc_dict
        self.out_dict.update(out_dict)
                    


if __name__ == '__main__':
    eft = EFTParam()
    eft.save_to_dict(year='2016', force_year='2018')
    eft.save_to_dict(year='2017', force_year='2018')
    eft.save_to_dict(year='2018', force_year='2018')
    import pickle
    #np.save(sys.path[1]+'Higgs-Combine-Tool/EFT_Parameterization_v1.npy', eft.out_dict, allow_pickle=True, fix_imports=True) # save to Higgs-Combine-Tool dir
    out_path = sys.path[1]+'Higgs-Combine-Tool/'
    out_file = 'EFT_Parameterization_v2.npy'
    #pickle.dump(eft.out_dict, open(out_path+out_file,'wb'), protocol=2) # save to Higgs-Combine-Tool dir
    # for ttbb testing
    out_path = sys.path[1]+'test/'
    out_file = 'EFT_Parameterization_ttbb.npy'
    pickle.dump(eft.eftbkg_fit, open(out_path+out_file,'wb'), protocol=2) # save to Higgs-Combine-Tool dir
    
    
