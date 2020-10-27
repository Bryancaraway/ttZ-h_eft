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
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import time
from collections import defaultdict
from numba import njit, prange
import concurrent.futures
import functools
#
import config.ana_cff as cfg
from lib.fun_library import fillne, t2Run
from modules.AnaDict import AnaDict
from modules.AnaVars import AnaVars
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
    outDir     = cfg.master_file_path
    njets      = 4
    maxAk4Jets = 100
    minBottoms = 2
    treeDir    = cfg.tree_dir
    jec_sys   = None
    jec_type  = None # ak4 or ak8 , might need to add sj for sdM
    doLepSel   = True
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
            #
            if 'UEUp' in self.sample and self.sample not in tree_dir: self.sample = self.sample.replace('UEUp','UEDUp') # because of some stupid type
            try:
                t = tree_dir.get(self.sample)
            except:
                print(f'{self.sample} not in input file for year: {self.year}\n Exitting...')
                exit()
            #
            self.DF_Container.set_tree(t, self.estart, self.estop)
            self.DF_Container.set_attr(self.isData, self.year, self.njets, self.maxAk4Jets, self.minBottoms, self.jec_sys, self.jec_type)
            if not self.isData: self.DF_Container.corr_AK8_jets() #un-disable this!
            self.DF_Container.set_current_tree_mask()
            if self.doLepSel: self.DF_Container.lep_sel_mask()
            # 
            ak4_df   = self.DF_Container('ak4',   'AK4_Variables')
            ak8_df   = self.DF_Container('ak8',   'AK8_Variables')
            val_df   = self.DF_Container('event', 'Event_Variables')
            ##
            val_df.df['n_ak4jets'], val_df.df['n_ak8jets'] = (self.DF_Container.n_ak4jets[self.DF_Container.event_mask], 
                                                              self.DF_Container.n_ak8jets[self.DF_Container.event_mask])
            if self.doLepSel:
                for k,v in self.DF_Container.lep_info.items():
                    val_df.df[k] = v 
            #
            rtc_df = self.cleanRC(ak4_df.df.loc(cfg.ana_vars['ak4lvec']['TLVarsLC']).sum())
            
            if self.isData:
                finish = time.perf_counter()
                print(f'\nTime to finish {type(self).__name__} for {self.sample}: {finish-start:.1f}\n')
                return ak4_df.df, ak8_df.df, val_df.df, rtc_df
            else:
                gen_df   = self.DF_Container('gen', 'GenPart_Variables')
                #TroubleShooting these systematics
                sys_df   = self.compute_ps_sc_weights().set_index(val_df.df.index)
                val_df.df = pd.concat([val_df.df, sys_df], axis='columns')
                # to-do add these to val_df
            #
            finish = time.perf_counter()
            print(f'\nTime to finish {type(self).__name__} for {self.sample}: {finish-start:.1f}\n')
            return ak4_df.df, ak8_df.df, val_df.df, rtc_df, gen_df.df
            #
    #
    @t2Run
    def cleanRC(self,ak4_LC):
        RC_lc_index = self.DF_Container('other', 'RC_LC_map', var=['Jet_lepcleaned_idx'])
        RC_vars   = self.DF_Container('RC',      'RC_TopInfo' )
        RC_vars   = RC_vars.df
        ##
        mapLC = fillne(RC_lc_index.df['Jet_lepcleaned_idx'])
        del RC_lc_index
        mapLC = np.where(mapLC==-1,np.nan,mapLC)
        RCj1, RCj2, RCj3 = fillne(RC_vars['ResolvedTopCandidate_j1Idx']), fillne(RC_vars['ResolvedTopCandidate_j2Idx']), fillne(RC_vars['ResolvedTopCandidate_j3Idx'])
        ##
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
    @t2Run
    def compute_ps_sc_weights(self):
        sys_df   = self.DF_Container('other', 'LHE_PS_weights', ['LHEScaleWeight', 'PSWeight','weight','LHEReweightingWeight'])
        #| Float_t LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; 
        #                                                           [3] is renscfact=1d0   facscfact=0.5d0 ; [4] is renscfact=1d0   facscfact=1d0 ; [5] is renscfact=1d0   facscfact=2d0 ; 
        #                                                           [6] is renscfact=2d0   facscfact=0.5d0 ; [7] is renscfact=2d0   facscfact=1d0 ; [8] is renscfact=2d0   facscfact=2d0 
        #| Float_t PS weights (w_var / w_nominal); [0] is ISR=0.5 FSR=1; [1] is ISR=1 FSR=0.5; [2] is ISR=2 FSR=1; [3] is ISR=1 FSR=2 *
        ps_w  = sys_df.df['PSWeight'].pad(4).fillna(1)
        sc_w  = sys_df.df['LHEScaleWeight'].pad(9).fillna(1)
        mc_w  = sys_df.df['weight']
        eft_w = sys_df.df['LHEReweightingWeight'].pad(184).fillna(1)
        #
        df = pd.DataFrame()
        # TROUBLESHOOT WHICH SAMPLES HAVE SCALE WEIGHT, PSWEIGHT

        df['ISR_Up']   = ps_w[:,2]
        df['ISR_Down'] = ps_w[:,0]
        df['FSR_Up']   = ps_w[:,3]
        df['FSR_Down'] = ps_w[:,1]
        #
        df['mu_r_Up']    = sc_w[:,7]
        df['mu_r_Down']  = sc_w[:,1]
        df['mu_f_Up']    = sc_w[:,5]
        df['mu_f_Down']  = sc_w[:,3]
        df['mu_rf_Up']   = sc_w[:,8]
        df['mu_rf_Down'] = sc_w[:,0]
        #
        for i in range(184):
            df[f'EFT{i}'] = eft_w[:,i]

        
        return df
        
    #
    class DF_Container: 
        '''
        container to dynamically handle root variables
        and save to .pkl files for further analysis
        container type must be ak4, ak8, event, gen, RC, or other
        '''
        tarray       = None
        tpandas      = None
        estart       = None
        estop        = None
        allowed_types = ['ak4', 'ak8', 'event', 'gen', 'RC', 'other']
    
        # pre-selection cuts need to be applied when getting data from tree to increase speed
        # object cuts: ak4jet_pt >= 30, |ak4jet_eta| <= 2.4 ... |ak8jet_eta| <= 2.4
        # event cuts:  n_ak4jets >= 4, n_ak4bottoms >= 2, n_ak8jets(pt>=200) >= 1
        isData     = False
        year       = None
        minAk4Jets = 4
        maxAk4Jets = 99
        minBottoms = 2
        #
        fj = {}
        #
        ak4_mask   = None
        ak8_mask   = None
        rtcd_mask  = None
        event_mask = None
        #
        n_ak4jets  = None
        n_ak8jets  = None
        
        def __init__(self, var_type, name, var=None,):
            #self.variables = (self.var_dict(var_type) if var is None else var)
            self.variables = (self.ana_vars[var_type] if var is None else var)
        
            self.var_type  = var_type
            self.name      = name
            # handle df depending on type 
            self.df = self.extract_and_cut()


        #@t2Run
        def extract_and_cut(self):
            type_indexer = defaultdict(
                None,{'ak4':  lambda v: self.build_dict(v)[self.ak4_mask][self.event_mask],
                      'ak8':  lambda v: AnaDict({**self.build_dict_ak8(v[:-2])[self.ak8_mask][self.event_mask],# super hacky but works
                                                 **self.build_dict(v[-2:])[self.event_mask]}), # fix in future, prone to break
                      'event':lambda v: self.tpandas(v)[self.event_mask],
                      'gen'  :lambda v: self.build_dict(v)[self.event_mask],
                      'RC'   :lambda v: self.build_dict(v)[self.rtcd_mask][self.event_mask],
                      'other':lambda v: self.build_dict(v)[self.event_mask]})
            try:
                df = type_indexer[self.var_type](self.variables)
            except KeyError:
                raise KeyError(f"At least one variable not found:{self.variables} \n Name '{self.var_type}' is not defined, Required to be: {self.allowed_types}")
            return df
            
        def apply_event_mask(self,df):
            return df[self.event_mask[df.index.get_level_values('entry')].values]
            
        @classmethod
        def set_tree(cls,tree, estart=None, estop=None):
            ''' Set tree array method and tree pandas method '''
            print(estart, estop)
            cls.tarray  = functools.partial(tree.array,      entrystart=estart, entrystop=estop)
            cls.tpandas = functools.partial(tree.pandas.df , entrystart=estart, entrystop=estop)
            
        @classmethod
        def set_attr(cls,isData, year, njets, maxAk4Jets, minBottoms, jec_sys=None, jec_type=None):
            cls.isData = isData
            cls.year   = year
            cls.minAk4Jets = njets
            cls.maxAk4Jets = maxAk4Jets
            cls.minBottoms = minBottoms
            cls.ana_vars = AnaVars(year,isData, jec_sys=jec_sys, jec_type=jec_type)

        @classmethod
        def set_current_tree_mask(cls):
            #
            jet_pt_eta    = cls.build_dict(cfg.ana_vars['ak4lvec']['TLVarsLC'][:2])
            fatjet_pt_eta = cls.build_dict_ak8(cfg.ana_vars['ak8lvec']['TLVarsLC'][:2])
            nbottoms     = cls.tarray(cls.ana_vars['nBottoms_drLeptonCleaned'])
            met_pt       = cls.tarray(cls.ana_vars['MET_pt'])
            
            rtcd = cls.tarray(cls.ana_vars[cfg.ana_vars['valRCvars'][0]])
            j_pt_key, j_eta_key   = cfg.ana_vars['ak4lvec']['TLVarsLC'][:2]
            fj_pt_key, fj_eta_key = cfg.ana_vars['ak8lvec']['TLVarsLC'][:2]
            
            #
            cls.ak4_mask  = ((jet_pt_eta[j_pt_key] >= 30) & (abs(jet_pt_eta[j_eta_key]) <= 2.4))
            cls.ak8_mask  = (abs(fatjet_pt_eta[fj_eta_key]) <= 2.4)
            cls.rtcd_mask = (rtcd >= 0.70) # trained with 75, might expand if I do re-train
            #
            n_ak4jets , n_ak8jets = jet_pt_eta[j_pt_key][cls.ak4_mask].counts, fatjet_pt_eta[fj_pt_key][(fatjet_pt_eta[fj_pt_key] >= cfg.ZHptcut) & (cls.ak8_mask)].counts
            del jet_pt_eta, fatjet_pt_eta
            #
            event_mask = ((n_ak4jets >= cls.minAk4Jets) & 
                          (n_ak4jets <= cls.maxAk4Jets) &
                          (n_ak8jets >= 1) &
                          (nbottoms >= cls.minBottoms) &
                          (met_pt > 20))
            #handle HEM Vetor for 2018 Data
            if cls.isData and cls.year == '2018':
                passHEMVeto = cls.tarray(cls.ana_vars['SAT_Pass_HEMVeto_DataOnly_drLeptonCleaned'])
                event_mask = (event_mask & (passHEMVeto == True))
            cls.n_ak4jets , cls.n_ak8jets = n_ak4jets , n_ak8jets
            cls.event_mask = event_mask
            #

        @classmethod
        def lep_sel_mask(cls):
            mu_info = cls.build_dict(cfg.lep_sel_vars['muon'])
            el_info = cls.build_dict(cfg.lep_sel_vars['electron'])
            #
            mu_mask =cfg.lep_sel['muon'](mu_info)  
            el_mask =cfg.lep_sel['electron'][cls.year](el_info)
            #
            lep_event_mask = ((mu_mask[mu_mask].counts + el_mask[el_mask].counts) == 1)
            cls.event_mask = cls.event_mask & lep_event_mask
            mu_info = mu_info[mu_mask][cls.event_mask]
            el_info = el_info[el_mask][cls.event_mask]
            
            #
            from awkward import JaggedArray as aj
            get_lep_info = (lambda k : aj.concatenate([mu_info[f'Muon_{k}'], el_info[f'Electron_{k}']],axis=1).flatten())
            single_mu = (mu_info['Muon_pt'].counts == 1)
            single_el = (el_info['Electron_pt'].counts == 1)
            
            cls.lep_info = AnaDict({'Lep_pt'  : get_lep_info('pt'),
                                    'Lep_eta' : get_lep_info('eta'),
                                    'Lep_phi' : get_lep_info('phi'),
                                    'Lep_mass': get_lep_info('mass'),
                                    'passSingleLepMu'   : single_mu,
                                    'passSingleLepElec' : single_el
                                })
        #
        @classmethod
        @t2Run
        def corr_AK8_jets(cls):
            ak8_vars = cls.build_dict(cfg.ak8_sys_vars, with_interp=False) #tries to get jec variation before its computed (have to initially use tarray)
            from modules.jmeAK8 import JMEAK8
            ak8_fj_transformer = JMEAK8(cls.year)
            #
            ak8_fj_transformer.prepare_transformer()
            cls.fj = ak8_fj_transformer.transform_AK8(ak8_vars)
            del ak8_fj_transformer
            #ak8_fj_transformer.transform_SDM(ak8_vars)
            #print(cls.fj.keys())
            
            
        @classmethod
        def build_dict(cls,keys, with_interp=True):
            executor = concurrent.futures.ThreadPoolExecutor()
            if with_interp:
                return AnaDict({k: cls.tarray(cls.ana_vars[k], executor=executor) for k in keys})
            else:
                return AnaDict({k: cls.tarray(k, executor=executor) for k in keys})

        @classmethod
        @t2Run
        def build_dict_ak8(cls,keys):
            executor = concurrent.futures.ThreadPoolExecutor()
            _dict = AnaDict({})
            for k in keys:
                #print(k)
                #if cls.ana_vars[k] in cls.fj:
                try:
                    _dict[k] = cls.fj[cls.ana_vars[k]]
                except: 
                    _dict[k] = cls.tarray(cls.ana_vars[k], executor=executor)
            return _dict
##

if __name__ == '__main__':
    #
    getData_cfg = {'files': ['MC_2017'], 'sample': 'DY', 'outDir': cfg.skim_ZHbb_dir,
                   'njets':cfg.ZHbbFitMinJets, 'maxAK4Jets':cfg.ZHbbFitMaxJets, 
                   'treeDir':cfg.tree_dir+'_bb', 'getGenData':True, 'getak8var':True,
                   'estop': None}
    #
    _ = getData(getData_cfg).getdata()
             
