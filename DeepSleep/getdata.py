########################
### Get data         ###
### for analysis     ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#
import uproot
import sys
import os
import math
import pickle
import operator
import re
import time
from collections import defaultdict
from numba import njit, prange
import concurrent.futures
import deepsleepcfg as cfg
#
import numpy as np
np.random.seed(0)
import pandas as pd
from itertools import combinations
##


class getData : 
    ''' 
    callable object designed to get relevent data from root file
    and transform into pandas dataframes 
    '''
    # Allowed class variables
    files      = cfg.files
    samples    = cfg.MCsamples
    outDir     = cfg.skim_dir
    njets      = 2
    maxAk4Jets  = 6
    treeDir    = cfg.tree_dir
    getGenData = False
    getak8var  = False
    
    

    def __init__ (self, kwargs):
        [setattr(self,k,v) for k,v in kwargs.items() if k in getData.__dict__.keys()]
        self.get_root_data()

    def get_root_data(self):
        # open root file
        for file_ in self.files: 
            if not os.path.exists(cfg.file_path+file_+'.root') : continue
            with uproot.open(cfg.file_path+file_+'.root') as f_:
                print(f'Opening File:\t{file_}')
                tree_dir = f_.get(self.treeDir)
                for sample in self.samples:
                    print(sample)
                    start = time.perf_counter()
                    t = tree_dir.get(sample)
                    self.DF_Container.set_current_tree_mask(t)
                    # 
                    idx = pd.IndexSlice
                    ak4_df   = self.DF_Container(cfg.ak4lvec['TLVarsLC']+cfg.ak4vars,'ak4',   'AK4_Variables')
                    ak8_df   = self.DF_Container(cfg.ak8lvec['TLVarsLC']+cfg.ak8vars,'ak8',   'AK8_Variables')
                    val_df   = self.DF_Container(cfg.valvars+cfg.sysvars,            'other', 'Event_Variables')
                    gen_df   = self.DF_Container(cfg.genpvars,                       'other', 'GenPart_Variables')
                    ##
                    val_df.df['n_ak4jets'], val_df.df['n_ak8jets'] = self.DF_Container.n_ak4jets, self.DF_Container.n_ak8jets, 
                    jet_df = pd.concat([ak4_df.df, ak8_df.df], keys=[ak4_df.name,ak8_df.name], names=['Jet_Type'])
                    #
                    rtc_df = self.cleanRC(jet_df.loc[ak4_df.name][cfg.ak4lvec['TLVarsLC']].sum(axis=1).unstack().to_numpy())
                    #
                    jet_df   .to_pickle('test/'+sample+file_[-5:]+'_jet.pkl')# replace test at some point FIXME
                    val_df.df.to_pickle('test/'+sample+file_[-5:]+'_val.pkl')# replace test at some point FIXME
                    gen_df.df.to_pickle('test/'+sample+file_[-5:]+'_gen.pkl')# replace test at some point FIXME
                    rtc_df   .to_pickle('test/'+sample+file_[-5:]+'_rtc.pkl')# replace test at some point FIXME
                    #
                    finish = time.perf_counter()
                    print(f'\nTime to finish {sample}: {finish-start:.1f}\n')
                #
    #
    def cleanRC(self,ak4_LC):
        RC_ak4    = self.DF_Container(cfg.ak4lvec['TLVars'], 'other',   'RC_AK4LVec' )
        RC_vars   = self.DF_Container(cfg.valRCvars,         'RC',   'RC_TopInfo' )
        RC_vars = RC_vars.df.unstack().reindex(self.DF_Container.event_mask[self.DF_Container.event_mask].index, fill_value=np.nan)
        ##
        ak4    = RC_ak4.df.sum(axis=1).unstack().to_numpy()
        ###
        RCj1, RCj2, RCj3 = RC_vars['ResolvedTopCandidate_j1Idx'].to_numpy(), RC_vars['ResolvedTopCandidate_j2Idx'].to_numpy(), RC_vars['ResolvedTopCandidate_j3Idx'].to_numpy()
        ##
        mapLC = self.cleanMap(ak4_LC, ak4, np.full(ak4.shape,np.nan))
        RCj1, RCj2, RCj3 = self.applyMap(mapLC, RCj1, RCj2, RCj3, np.full(RCj1.shape,np.nan))
        RC_vars['ResolvedTopCandidate_j1Idx'], RC_vars['ResolvedTopCandidate_j2Idx'], RC_vars['ResolvedTopCandidate_j3Idx'] = RCj1, RCj2, RCj3
        return RC_vars.stack(dropna=False)

    @staticmethod
    @njit(parallel=True)
    def cleanMap(a,b, out):
        rows, a_cols = a.shape
        b_cols = b.shape[1]
        for i in prange(rows):
            for j in range(a_cols):
                if np.isnan(a[i,j]):
                    continue
                for k in range(b_cols):
                    match = a[i,j] == b[i,k]
                    if match:
                        out[i,k] = int(j)
                        break
        return out

    @staticmethod
    @njit(parallel=True)
    def applyMap(lc, j1, j2, j3, o1):
        rows, n_jcols = j1.shape
        o2 = o3 = o1
        for i in prange(rows):
            for j in range(n_jcols):
                if np.isnan(j1[i,j]):
                    continue
                id1 = int(j1[i,j])
                id2 = int(j2[i,j])
                id3 = int(j3[i,j])
                o1[i,j] = lc[i,id1]
                o2[i,j] = lc[i,id2]
                o3[i,j] = lc[i,id3]
        return o1, o2, o3
    #
    class DF_Container: 
        '''
        container to dynamically handle root variables
        and save to .pkl files for further analysis
        container type must be ak4, ak8, or other
        '''
        current_tree = None
        allowed_types = ['ak4', 'ak8', 'other', 'RC']
    
        # pre-selection cuts need to be applied when getting data from tree to increase speed
        # object cuts: ak4jet_pt >= 30, |ak4jet_eta| <= 2.4 ... |ak8jet_eta| <= 2.0
        # event cuts:  n_ak4jets >= 4, n_ak4bottoms >= 2, n_ak8jets(pt>=200) >= 1
        ak4_mask   = None
        ak8_mask   = None
        rtcd_mask  = None
        event_mask = None
        #
        n_ak4jets  = None
        n_ak8jets  = None
        
        def __init__(self, variables, var_type, name):
            self.variables = variables
            self.var_type  = var_type
            self.name      = name
            # handle df depending on type 
            self.df = self.extract_and_cut()
    
        def extract_and_cut(self):
            idx = pd.IndexSlice
            tree_to_df = self.current_tree.pandas.df
            type_indexer = defaultdict(
                None,{'ak4':  lambda v: self.apply_event_mask(tree_to_df(v)[self.ak4_mask]),
                      'ak8':  lambda v: pd.concat([self.apply_event_mask(tree_to_df(v[:-2])[self.ak8_mask]),  # seperate subjet in config FIXME
                                                  self.apply_event_mask(tree_to_df(v[-2:]))], axis='columns'),# seperate subjet in config FIXME 
                      'other':lambda v: self.apply_event_mask(tree_to_df(v)),
                      'RC':   lambda v: self.apply_event_mask(tree_to_df(v)[self.rtcd_mask])})
            try:
                df = type_indexer[self.var_type](self.variables)
            except KeyError:
                raise KeyError(f"Name '{self.var_type}' is not defined, Required to be: {self.allowed_types}")
            return df
            
        def apply_event_mask(self,df):
            return df[self.event_mask[df.index.get_level_values('entry')].values]
    
        @classmethod
        def set_current_tree_mask(cls,tree):
            cls.current_tree = tree
            #
            jet_pt_eta    = tree.pandas.df(cfg.ak4lvec['TLVarsLC'][:2])
            fatjet_pt_eta = tree.pandas.df(cfg.ak8lvec['TLVarsLC'][:2])
            rtcd = tree.pandas.df(cfg.valRCvars[0])
            jet_pt_key, jet_eta_key = list(jet_pt_eta.columns)
            fatjet_pt_key, fatjet_eta_key = list(fatjet_pt_eta.columns)
            #
            cls.ak4_mask  = ((jet_pt_eta[jet_pt_key] >= 30) & (abs(jet_pt_eta[jet_eta_key]) <= 2.4))
            cls.ak8_mask  = (abs(fatjet_pt_eta[fatjet_eta_key]) <= 2.0)
            cls.rtcd_mask = (rtcd[cfg.valRCvars[0]] >= 0.75)
            #
            cls.n_ak4jets , cls.n_ak8jets = jet_pt_eta[jet_pt_key][cls.ak4_mask].count(level='entry'), fatjet_pt_eta[fatjet_pt_eta>=200][fatjet_pt_key][cls.ak8_mask].count(level='entry')
            del jet_pt_eta, fatjet_pt_eta
            #
            cls.event_mask = ((cls.n_ak4jets >= getData.njets) & 
                              (cls.n_ak4jets <= getData.maxAk4Jets) &
                              (cls.n_ak8jets >= 1))
            
##

if __name__ == '__main__':
    #
    getData_cfg = {'files': ['result_2017'], 'samples': cfg.ZHbbFitCfg[1], 'outDir': cfg.skim_ZHbb_dir,
                   'njets':cfg.ZHbbFitCut[1], 'maxJets':cfg.ZHbbFitMaxJets, 
                   'treeDir':cfg.tree_dir+'_bb', 'getGenData':True, 'getak8var':True}
    #
    getData(getData_cfg)

             
