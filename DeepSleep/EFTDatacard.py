
#                      #  
##                    ##
########################                               
### TTZ/H, Z/H to bb ###
### build datacard   ###                               
### for EFT, HCT     ###                               
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
#import pickle
#import math
import functools
#from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Pool
import re
#import operator as OP
#
import config.ana_cff as cfg
#import config.process_norms as p_norms
from lib.fun_library import weighted_quantile, getZhbbBaseCuts, getZhbbWeight, t2Run
from lib.TH1 import export1d
#from lib.datacard_shapes import DataCardShapes
from modules.eftParam import EFTParam
#
import uproot
#import uproot_methods
#from ROOT import TFile, TDirectory, TH1F

import numpy as np
#import pandas as pd
#import json
#import matplotlib.pyplot as plt

class EFTDatacard:
    '''
    This class grabs existing datacard
    and converts the shapes to channels
    '''
    tag = '' if len(sys.argv) < 2 else sys.argv[1]+'_'
    tag = 'recoeft_'
    dc_dir = 'Higgs-Combine-Tool'
    eft_dir = 'Higgs-Combine-Tool/eftdatacards'
    dummy_hist = 'ch_TTBar'
    init_ch = ['Zhpt1','Zhpt2','Zhpt3']
    #eft_out_file = 'EFT_Parameterization_v4.npy'
    def __init__(self):
        self.years = cfg.Years
        self.datacards = [f'{self.dc_dir}/datacard_{self.tag}{y}.txt' for y in self.years]
        self.dataroos  = [f'{self.dc_dir}/datacard_{self.tag}{y}.root' for y in self.years]
        self.ch_bins = {
            ch: len(uproot.open(self.dataroos[0])[self.dummy_hist.replace('ch',ch)].values) for ch in self.init_ch
        }
        self.max_ch_bins  = max([self.ch_bins[chb] for chb in self.ch_bins])
        self.n_ch = len(self.init_ch)
        #
        # make new datacards and root file
        self.copy_create_ch_txt() 
        self.copy_create_ch_root()
        
    @t2Run
    def copy_create_ch_txt(self):
        for datacard in self.datacards:
            with open(datacard,'r') as r_txt:
                old_dc = r_txt.readlines()
                for i_bin in range(self.max_ch_bins):
                    new_dc_name = re.sub(r'(201\d)', rf'{i_bin}_\1',datacard.replace(self.dc_dir, self.eft_dir))
                    with open(new_dc_name,'w') as out_txt:
                        if i_bin < self.ch_bins['Zhpt1']:
                            out_txt.writelines([re.sub(r'(Zhpt\d)',rf'\1_{i_bin}',l) for l in old_dc])
                        else:
                            self.handle_non_uniform_dc(old_dc,out_txt,i_bin)


    def handle_non_uniform_dc(self, old, new, i_bin):
        _out_txt = []
        for line in old:
            offset = (1 if 'lnN' not in line and 'shape' not in line else 2)
            n_proc = (len(line.split())-offset)
            if n_proc % self.n_ch == 0  and 'jmax' not in line and '---' not in line and 'shapes' not in line:
                l_index = int(n_proc/self.n_ch)
                pre_line = line.split()[:offset]
                out_line = '  '.join(pre_line+line.split()[l_index+offset:])+' \n'
                _out_txt.append(re.sub(r'(Zhpt\d)', rf'\1_{i_bin}', out_line))
            else:
                _out_txt.append(line)
        #
        new.writelines(_out_txt)
        

    @t2Run
    def copy_create_ch_root(self):
        pool = Pool()
        _ = pool.map(self.roo_worker, self.dataroos)
        pool.close()
        
    def roo_worker(self,dataroo):
        with uproot.open(dataroo) as roo:
            roo_name = dataroo.replace(self.dc_dir,self.eft_dir)
            if os.path.exists(roo_name):
                os.system(f"rm {roo_name}")
            new_roo = uproot.create(roo_name)
            _out_dict = {}
            for hist in roo:
                for i,(val,var) in enumerate(zip(roo[hist].values,roo[hist].variances)):
                    #if (val == 0 or var == 0) and 'Up' not in hist.decode() and 'Down' not in hist.decode():
                    #if (val == 0 ) :
                    #    print(f'Zero bin entry in {roo} for {hist}, bin {i} : {val},{var}')
                    _dict = {'sumw':np.array([val]),'sumw2':np.array([var])}
                    _name = re.sub(r';\d*','',re.sub(r'(Zhpt\d)',rf'\1_{i}',hist.decode()))
                    _out_dict[_name] = export1d(_dict,_name, z_to_e=True)
            #
            for k in _out_dict:
                new_roo[k] = _out_dict[k]
            new_roo.close()
                    
def write2roo(k,outdict=None,outroo=None):
    print(k)
    outroo[k] = outdict[k]

if __name__ == '__main__':
    EFTDatacard()
