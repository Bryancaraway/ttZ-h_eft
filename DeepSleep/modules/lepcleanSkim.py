########################
### Skim data        ###
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
import json
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
from lib.fun_library import fillne, t2Run, deltaR, ak_crosscleaned
from modules.AnaDict import AnaDict
from modules.AnaVars import AnaVars
from modules.Skim import Skim
from modules.metaSkim import SkimMeta
#
import numpy as np
from awkward import JaggedArray as aj
np.random.seed(0)
import pandas as pd
##

class LepCleanSkim(Skim) :
    '''
    messed up skim so need to 
    fix current skim
    ----
    open skimmed pkl --> clean jets from leptons
    --> reapply cuts --> update metaData -- rewrite to pkl
    '''
    @t2Run
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        

    @t2Run
    def startSkim(self):
        self.open_pickle(self.roofile) # not a root file, but pickle
        del self.roofile
        self.Meta = SkimMeta(self.sample, self.year, self.isData, 
                             tree=None, jec_sys=self.jec_sys, run_norm=False)
        # define event information
        # apply in house jmr 
        self.fatjets   = self.build_dict('ak8') # jus in case jec
        #self.jets      = self.build_dict('ak4')
        self.events    = self.build_dict('events')
        self.jets = self.build_ak4jets()
        exit()
        self.electrons = AnaDict({
            k.replace('Lep','Electron'): aj.fromcounts(
                self.events['passSingleLepElec'].to_numpy(int),self.events[k]) 
            for k in ['Lep_pt','Lep_eta','Lep_phi']})
        self.muons     = AnaDict({
            k.replace('Lep','Muon'): aj.fromcounts(
                self.events['passSingleLepMu'].to_numpy(int),self.events[k]) 
            for k in ['Lep_pt','Lep_eta','Lep_phi']}) 
        # other things like gen info
        if not self.isData:
            self.geninfo    = self.build_dict('gen') 
        #self.subjets    = self.build_dict(cfg.ana_vars['ak8sj'])
        # wont keep
        # ===================== #
        # apply object criteria
        self.jet_dr_veto = self.is_a_jet()
        self.jets           = self.jets[     self.is_a_jet()] # still need to lepton clean
        self.fatjets        = self.fatjets[  self.is_a_fatjet()]
        # define event selection
        self.event_mask   = self.get_event_selection()
        # apply event selection
        self.jet_dr_veto     = self.jet_dr_veto[   self.event_mask]
        self.jets            = self.jets[          self.event_mask]
        self.fatjets         = self.fatjets[       self.event_mask]
        self.electrons       = self.electrons[     self.event_mask]
        self.muons           = self.muons[         self.event_mask]
        self.events          = self.events[        self.event_mask]
        if not self.isData:
            self.geninfo    = self.geninfo[   self.event_mask]
        ''' 
        add interesting info to events: 
        (lepton pt, eta, phi , etc...) 
        njets, nfjets, nbjets
        ps, sc, eft
        hem veto weight
        '''
        self.handle_multiplicity_HEM_info()
        #if not self.isData:
        #    need to open all files and get jet btag sf shapes
        #    self.Meta.add_btagweightsf_counts(self.jets, self.events, self.geninfo) # self.jet_dr_veto
        #
        self.btag_event_mask = self.get_b_tag_eventmask() 
        self.jets            = self.jets[          self.btag_event_mask]
        self.fatjets         = self.fatjets[       self.btag_event_mask]
        #self.subjets         = self.subjets[       self.btag_event_mask]
        self.events          = self.events [       self.btag_event_mask]
        if not self.isData:
            self.geninfo = self.geninfo[  self.btag_event_mask]
        #
        
    
    def get_updated_skim(self):
        __out_dict = self.f
        __out_dict['ak4'] = self.jets
        __out_dict['ak8'] = self.fatjets
        __out_dict['events'] = self.events
        if not self.isData:
            __out_dict['gen'] = self.geninfo
            __out_dict['metaData'] = self.Meta.get_metadata()
        return __out_dict

    def get_MET_filter(self):
        return True

    def open_pickle(self, ifile):
        ''' opening from pkl, not root '''
        self.f = AnaDict.read_pickle(ifile)
        
    def build_dict(self,  ktype, with_interp=True):
        # pass type: metaData, events, ak4, ak8
        #if kype != 'events':
        return self.f[ktype]
                
    @t2Run
    def build_ak4jets(self):
        '''
        need to open every file and 
        retreive events that match run, lumiblock, event #
        '''
        from glob import glob
        # does not work on data
        self.ana_vars = AnaVars(self.year, self.isData, jec_sys=self.jec_sys)
        executor = concurrent.futures.ThreadPoolExecutor()
        __ak4_dict = AnaDict({})
        rfiles = files = glob(f"{cfg.postproc_dir}/{self.year}/{self.sample}_{self.year}/*.root")
        with uproot.open(rfiles[0]) as roo:
            _tree = roo.get('Events')
            print(_tree.array(['nJet']))
            _rle       = _tree.pandas.df(['run','luminosityBlock','event'])
            _srle      = self.events.loc[:,['run','luminosityBlock','event']]
            _ind = _rle['run'].isin(_srle['run']) & \
                   _rle['luminosityBlock'].isin(_srle['luminosityBlock']) & \
                   _rle['event'].isin(_srle['event']) 
            print(_ind)
            _jets = AnaDict({
                k: _tree.array(self.ana_vars[k], executor=executor)[_ind] for k in
                cfg.ana_vars['ak4vars']+cfg.ana_vars['ak4lvec']['TLVars']+(
                    [] if self.isData else cfg.ana_vars['ak4mcvars']
                )})
            print(_jets)
        exit()
            
            
        

if __name__ == '__main__':
    # test quick lep clean
    from config.sample_cff import sample_cfg, process_cfg
    #sample = 'TTToSemiLeptonic'
    sample  = 'TTZToQQ'
    year = '2016'
    golden_json=json.load(open(cfg.goodLumis_file[year]))
    test_file = cfg.postSkim_dir+\
            f"{year}/{sample_cfg[sample]['out_name']}/{sample}.jesHFUp.pkl"
    
    _ = LepCleanSkim(test_file, sample, year, isData='Data' in sample, jec_sys='jesHFUp', golden_json=golden_json)
    out_file = _.get_updated_skim()
    start_file = AnaDict.read_pickle(test_file) 
    #print(sum(out_file['ak4']['Jet_pt'].counts == start_file['ak4']['Jet_pt'].counts))
    print(len(out_file['ak4']['Jet_pt'].counts))
    print(len(start_file['ak4']['Jet_pt'].counts))
