import os
import time
import argparse
import re
import numpy as np
import pandas as pd
import json
from glob import glob
import subprocess as sb
#
import config.ana_cff as cfg
from config.sample_cff import sample_cfg
from lib.fun_library import t2Run
#  module classes
from modules.AnaDict import AnaDict


class PostSkim :
    '''
    Script run at the end
    of Skimming
    '''
    def __init__(self, sample, year, isData, out_dir, tag):
        self.isData  = isData
        self.lumi    = cfg.Lumi[year]
        self.outfile = f"{out_dir}/{sample if not isData else sample_cfg[sample]['out_name']}{tag}.pkl"
        self.files   = (lambda : glob(f'{out_dir}/{sample}_*{tag}*.pkl'))
        print(out_dir)
        self.metaData  = {'sample':sample,'year':year,
                          'xs': sample_cfg[sample]['xs'],
                          'kf': sample_cfg[sample]['kf']}
        self.interp_dict = {
            'events'   : self.concat_events,
            'metaData' : self.concat_meta,
            # concat_other
        }
        self.final_pkl = {}

    def run(self):
        self.concat_files()
        #self.final_pkl['metaData'].update(self.metaData)
        print(self.final_pkl.keys())
        self.final_pkl['events']['sample'] = self.metaData['sample']
        if self.isData:
            weight = 1
        else:
            weight = (self.metaData['xs']*self.metaData['kf']*self.lumi*1000)/self.final_pkl['metaData']['tot_events']
            self.handle_btag_weight()
        self.final_pkl['events']['weight'] = weight
        AnaDict(self.final_pkl).to_pickle(self.outfile)

    # --- #
    def concat_files(self):
        self.files = self.files()
        self.final_pkl = self.open_and_del(self.files[0])
        for pkl in self.files[1:]:        
            pkl_dict = self.open_and_del(pkl)
            for k in pkl_dict:
                #print(pkl_dict[k])
                self.interp_dict.get(k, (lambda awk: self.concat_other(awk,k)))(pkl_dict[k])
    
    #
    def concat_events(self,df):
        self.final_pkl['events'] = pd.concat([self.final_pkl['events'],df], 
                                             axis = 'rows', ignore_index=True)
    def concat_meta(self,di):
        self.final_pkl['metaData'] = {var: self.final_pkl['metaData'][var] + di[var] for var in di}
    def concat_other(self,awk,k):
        self.final_pkl[k] = {var: self.try_concatenate([self.final_pkl[k][var],awk[var]]) for var in awk}
    # == #
    def handle_btag_weight(self):
        df = self.final_pkl['events'].filter(like='BTagWeight', axis='columns')
        nj = self.final_pkl['events']['n_ak4jets'].clip(0,12)
        for k in df.keys():
            bty = 'btw_yield'+k.replace('BTagWeight','')
            opp = np.where(self.final_pkl['metaData'] == 0 , 0.,
                           self.final_pkl['metaData']['nj_yield']/self.final_pkl['metaData'][bty])
            self.final_pkl['events'].loc[:,k] = df[k]*nj.apply((lambda i : opp[i]))
    # -- #
    @staticmethod
    def open_and_del(pkl):
        pkl_dict = AnaDict.read_pickle(pkl)
        os.system(f'rm {pkl}')
        return pkl_dict
    @staticmethod
    def try_concatenate(arrays):
        import awkward
        try:
            contents = np.concatenate([j.flatten() for j in arrays])
            counts = np.concatenate([j.counts for j in arrays])
            return awkward.JaggedArray.fromcounts(counts, contents)
        except AttributeError:
            return np.concatenate([arr for arr in arrays])
            


if __name__ == '__main__':
    PostSkim()
