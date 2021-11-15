##### Calculate R Ratio for BTagW #####
### Written by: Bryan Caraway       ###
#######################################

import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import numpy as np
import pandas as pd
import json
import re
#import functools
import config.ana_cff as cfg
from config.sample_cff import sample_cfg, process_cfg
from modules.AnaDict import AnaDict

class Calc_BTagWeightSF : 
    '''
    open metadata from skim
    and calculate r ratio by process
    '''
    
    SkimDir = f"{cfg.postSkim_dir}/"

    def __init__(self, year) : 
        self.year       = year
        self.processes  = cfg.All_MC+['QCD']+['VV']
        #self.btw_yields = None
        self.r_ratio    = {} # by process
        #
        self.concat_btw_meta()
        self.write_to_json()
    
    def concat_btw_meta(self):
        for p in self.processes:
            bty = None
            if p not in process_cfg: continue
            for mc in process_cfg[p]:
                mD = AnaDict.read_pickle(f'{self.SkimDir}/{self.year}/{p}/{mc}.pkl')['metaData']
                #if self.btw_yields is None: self.btw_yields = re.findall(r'\w*_yield\w*', ' '.join(mD.keys()))
                bty_names = re.findall(r'\w*_yield\w*', ' '.join(mD.keys()))
                bty = self.update_dict(mD,bty,bty_names) # in, out, keys
            #
            #if p == 'TTBar': # do tt+cc and tt+lf seperately
            #    self.r_ratio[f'{p}_lf'] = self.calc_r_ratio(bty, opt='lf')
            #    self.r_ratio[f'{p}_cc'] = self.calc_r_ratio(bty, opt='cc')
            #else:
            self.r_ratio[p] = self.calc_r_ratio(bty)

    def write_to_json(self):
        out_name = cfg.dataDir+f'/btagw_r_ratio/btagw_r_ratio_{self.year}.json'
        with open(out_name,'w') as jsf:
            json.dump(self.r_ratio, jsf, indent=4)

    @staticmethod
    def update_dict(in_,out_,keys_):
        if out_ is None: 
            return {k_ : in_[k_]*in_['weight'] for k_ in keys_}
        return {k_ : in_[k_]*in_['weight']+out_[k_] for k_ in keys_}
    
    @staticmethod
    def calc_r_ratio(yields_, opt=None):
        # clip at 9 
        opt      = f'{opt}_' if opt is not None else ''
        arr_clip =  (lambda arr,i: np.append(arr[:i],sum(arr[i:]))) 
        nj_yield = yields_.pop(f'{opt}nj_yield')
        nj_yield = arr_clip(nj_yield, 9)
        r_ratio = {}
        for y_ in yields_:
            if opt not in y_: continue
            bty_ = arr_clip(yields_[y_],9)
            r_ratio[y_.replace(f'{opt}btw_yield','r_ratio')] = np.where(bty_ == 0 , 0., nj_yield/bty_).tolist()
        #print(r_ratio)
        return r_ratio

if __name__ == '__main__':
    #
    Calc_BTagWeightSF('2016')
    Calc_BTagWeightSF('2017')
    Calc_BTagWeightSF('2018')
