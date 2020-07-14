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
import os
import re
import time
from collections import defaultdict
from numba import njit, prange
import concurrent.futures
import functools
#
from cfg import deepsleepcfg as cfg
from lib.fun_library import fillne, t2Run
from modules.AnaDict import AnaDict
#
import numpy as np
np.random.seed(0)
import pandas as pd
##


class getData : 
    ''' 
    callable object designed to get relevent data from root file
    and transform into pandas dataframes 
    '''
    # Allowed class variables
    isData     = False
    roofile    = None
    year       = None
    sample     = None
    outDir     = 'files/'
    njets      = 4
    maxAk4Jets = 100
    treeDir    = cfg.tree_dir
    getGenData = False # doesnt do anything atm
    getak8var  = False # doesnt do anything atm
    estart     = None
    estop      = None
    

    def __init__ (self, kwargs=None):
        [setattr(self,k,v) for k,v in kwargs.items() if k in getData.__dict__.keys()]

    def getdata(self):
        # open root file
        #if not os.path.exists(cfg.file_path+self.roofile+'.root') : raise
        with uproot.open(cfg.file_path+self.roofile+'.root') as f_:
            print(f'Opening File:\t{self.roofile}')
            tree_dir = f_.get(self.treeDir)
            print(self.sample)
            start = time.perf_counter()
            t = tree_dir.get(self.sample)
            self.DF_Container.set_attr(self.isData, self.year, self.njets, self.maxAk4Jets, self.estart, self.estop)
            self.DF_Container.set_current_tree_mask(t)
            # 
            ak4_df   = self.DF_Container('ak4',   'AK4_Variables')
            ak8_df   = self.DF_Container('ak8',   'AK8_Variables')
            val_df   = self.DF_Container('event', 'Event_Variables')
            ##
            val_df.df['n_ak4jets'], val_df.df['n_ak8jets'] = self.DF_Container.n_ak4jets, self.DF_Container.n_ak8jets
            #
            rtc_df = self.cleanRC(ak4_df.df.loc(cfg.ana_vars['ak4lvec']['TLVarsLC']).sum())
            #
            if self.isData:
                finish = time.perf_counter()
                print(f'\nTime to finish {type(self).__name__} for {self.sample}: {finish-start:.1f}\n')
                return ak4_df.df, ak8_df.df, val_df.df, rtc_df
            else:
                gen_df   = self.DF_Container('gen', 'GenPart_Variables')
                finish = time.perf_counter()
                print(f'\nTime to finish {type(self).__name__} for {self.sample}: {finish-start:.1f}\n')
                return ak4_df.df, ak8_df.df, val_df.df, rtc_df, gen_df.df
            #
    #
    #@t2Run
    def cleanRC(self,ak4_LC):
        RC_ak4    = self.DF_Container('other',   'RC_AK4LVec', var=cfg.ana_vars['ak4lvec']['TLVars'])
        RC_vars   = self.DF_Container('RC',      'RC_TopInfo' )
        RC_vars   = RC_vars.df
        ##
        ak4    = fillne(RC_ak4.df.sum())
        RCj1, RCj2, RCj3 = fillne(RC_vars['ResolvedTopCandidate_j1Idx']), fillne(RC_vars['ResolvedTopCandidate_j2Idx']), fillne(RC_vars['ResolvedTopCandidate_j3Idx'])
        ##
        mapLC = self.cleanMap(fillne(ak4_LC), ak4, np.full(ak4.shape,np.nan))
        RCj1, RCj2, RCj3 = self.applyMap(mapLC, RCj1, RCj2, RCj3, np.full(RCj1.shape,np.nan), np.full(RCj1.shape,np.nan), np.full(RCj1.shape,np.nan))
        RC_vars['ResolvedTopCandidate_j1Idx'], RC_vars['ResolvedTopCandidate_j2Idx'], RC_vars['ResolvedTopCandidate_j3Idx'] = RCj1, RCj2, RCj3
        return RC_vars

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
    #

    @staticmethod
    @njit(parallel=True)
    def applyMap(lc, j1, j2, j3, o1, o2, o3):
        rows, n_jcols = j1.shape
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
        container type must be ak4, ak8, event, gen, RC, or other
        '''
        tarray       = None
        tpandas      = None
        estop        = None
        allowed_types = ['ak4', 'ak8', 'event', 'gen', 'RC', 'other']
    
        # pre-selection cuts need to be applied when getting data from tree to increase speed
        # object cuts: ak4jet_pt >= 30, |ak4jet_eta| <= 2.4 ... |ak8jet_eta| <= 2.4
        # event cuts:  n_ak4jets >= 4, n_ak4bottoms >= 2, n_ak8jets(pt>=200) >= 1
        isData     = False
        year       = None
        minAk4Jets = 4
        maxAk4Jets = 99

        ak4_mask   = None
        ak8_mask   = None
        rtcd_mask  = None
        event_mask = None
        #
        n_ak4jets  = None
        n_ak8jets  = None
        
        def __init__(self, var_type, name, var=None,):
            self.variables = (self.var_dict(var_type) if var is None else var)
        
            self.var_type  = var_type
            self.name      = name
            # handle df depending on type 
            self.df = self.extract_and_cut()

            
        def var_dict(self,_type):
            _dict  = {
                'ak4'   : cfg.ana_vars['ak4lvec']['TLVarsLC']+cfg.ana_vars['ak4vars'],
                'ak8'   : cfg.ana_vars['ak8lvec']['TLVarsLC']+cfg.ana_vars['ak8vars']+cfg.ana_vars['ak8sj'],
                'event' : cfg.ana_vars['valvars']+(['run']+cfg.ana_vars['dataHLT_all']+cfg.ana_vars[f'dataHLT_{self.year}'] if self.isData else 
                                                   cfg.ana_vars['sysvars_mc']+cfg.ana_vars[f'sysvars_{self.year}'])+
                (cfg.ana_vars['HEM_veto'] if self.year == '2018' else []),
                'gen'   : cfg.ana_vars['genpvars'],
                'RC'    : cfg.ana_vars['valRCvars']}
            return _dict[_type]

        #@t2Run
        def extract_and_cut(self):
            type_indexer = defaultdict(
                None,{'ak4':  lambda v: self.build_dict(v)[self.ak4_mask][self.event_mask],
                      'ak8':  lambda v: AnaDict({**self.build_dict(v[:-2])[self.ak8_mask][self.event_mask],# super hacky but works
                                                 **self.build_dict(v[-2:])[self.event_mask]}), # fix in future, prone to break
                      'event':lambda v: self.tpandas(v)[self.event_mask],
                      'gen'  :lambda v: self.build_dict(v)[self.event_mask],
                      'RC'   :lambda v: self.build_dict(v)[self.rtcd_mask][self.event_mask],
                      'other':lambda v: self.build_dict(v)[self.event_mask]})
            try:
                df = type_indexer[self.var_type](self.variables)
            except KeyError:
                raise KeyError(f"At least one variable not found:{self.variables}\nName '{self.var_type}' is not defined, Required to be: {self.allowed_types}")
            return df
            
        def apply_event_mask(self,df):
            return df[self.event_mask[df.index.get_level_values('entry')].values]
            
        @classmethod
        def set_attr(cls,isData, year, njets, maxAk4Jets, estart, estop):
            cls.isData = isData
            cls.year   = year
            cls.minAk4Jets = njets
            cls.maxAk4Jets = maxAk4Jets
            cls.estart = estart
            cls.estop  = estop

        @classmethod
        def set_current_tree_mask(cls,tree):
            cls.tarray  = functools.partial(tree.array,      entrystop=cls.estop)
            cls.tpandas = functools.partial(tree.pandas.df , entrystop=cls.estop)
            #
            jet_pt_eta    = cls.build_dict(cfg.ana_vars['ak4lvec']['TLVarsLC'][:2])
            fatjet_pt_eta = cls.build_dict(cfg.ana_vars['ak8lvec']['TLVarsLC'][:2])
            
            rtcd = cls.tarray(cfg.ana_vars['valRCvars'][0])
            j_pt_key, j_eta_key   = cfg.ana_vars['ak4lvec']['TLVarsLC'][:2]
            fj_pt_key, fj_eta_key = cfg.ana_vars['ak8lvec']['TLVarsLC'][:2]
            
            #
            cls.ak4_mask  = ((jet_pt_eta[j_pt_key] >= 30) & (abs(jet_pt_eta[j_eta_key]) <= 2.4))
            cls.ak8_mask  = (abs(fatjet_pt_eta[fj_eta_key]) <= 2.4)
            cls.rtcd_mask = (rtcd >= 0.50)
            #
            n_ak4jets , n_ak8jets = jet_pt_eta[j_pt_key][cls.ak4_mask].counts, fatjet_pt_eta[fj_pt_key][(fatjet_pt_eta[fj_pt_key] >= cfg.ZHptcut) & (cls.ak8_mask)].counts
            del jet_pt_eta, fatjet_pt_eta
            #
            event_mask = ((n_ak4jets >= cls.minAk4Jets) & 
                          (n_ak4jets <= cls.maxAk4Jets) &
                          (n_ak8jets >= 1))
            #handle HEM Vetor for 2018 Data
            if cls.isData and cls.year == '2018':
                passHEMVeto = cls.tarray('SAT_Pass_HEMVeto_DataOnly_drLeptonCleaned')
                event_mask = (event_mask & (passHEMVeto == True))
            cls.n_ak4jets , cls.n_ak8jets = n_ak4jets[event_mask] , n_ak8jets[event_mask]
            cls.event_mask = event_mask
            #
        @classmethod
        def build_dict(cls,keys):
            executor = concurrent.futures.ThreadPoolExecutor()
            return AnaDict({k: cls.tarray(k, executor=executor) for k in keys})
            

##

if __name__ == '__main__':
    #
    getData_cfg = {'files': ['MC_2017'], 'sample': 'DY', 'outDir': cfg.skim_ZHbb_dir,
                   'njets':cfg.ZHbbFitMinJets, 'maxAK4Jets':cfg.ZHbbFitMaxJets, 
                   'treeDir':cfg.tree_dir+'_bb', 'getGenData':True, 'getak8var':True,
                   'estop': None}
    #
    _ = getData(getData_cfg).getdata()
             
