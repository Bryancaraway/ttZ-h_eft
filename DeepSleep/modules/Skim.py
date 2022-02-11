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
from modules.metaSkim import SkimMeta
from modules.ak8jmsjmr_helper import ak8jmsjmr_helper
from modules.pdfweight_helper import PDFHelper
#
import numpy as np
from awkward import JaggedArray as aj
np.random.seed(0)
import pandas as pd
##

class Skim :
    '''
    This code does two things:
    1. define analysis objects: leptons, jets, fatjets, bjets
    2. apply event level cuts/filters to elliminate unwanted events
    3. (for data) apply golden json file
    '''
    @t2Run
    def __init__(self, roofile, sample, year, isData=False, is4eff=False, jec_sys=None, golden_json=None):
        self.roofile = roofile
        print(roofile) # for potential debugging 
        self.sample = sample
        self.year   = year
        self.isData = isData
        self.is4eff = is4eff
        self.jec_sys = jec_sys
        #
        self.golden_json = golden_json
        #
        self.startSkim()
    @t2Run
    def startSkim(self):
        self.ana_vars = AnaVars(self.year, self.isData, jec_sys=self.jec_sys) 
        self.tree     = self.set_tree_from_roofile(self.roofile)
        # handle pdfweight calc/matching if needed
        #if not self.isData:
        pdf_helper = PDFHelper(self) # currently does nothing for EFT samples
        # prepMeta metadata factory
        self.Meta = SkimMeta(self.sample, self.year, self.isData, self.tree, self.jec_sys, pdf=pdf_helper.pdfweights)
        # define event information
        # apply in house jmr 
        self.fatjets   = self.build_dict(['FatJet_pt','FatJet_msoftdrop']) # just in case jec
        self.tmp_fatjets = self.build_dict(cfg.ana_vars['ak8vars']+cfg.ana_vars['ak8lvec']['TLVars_nom'], with_interp=False)
        self.tmp_fatjets['FatJet_pt'] = self.tmp_fatjets.pop('FatJet_pt_nom') # rename this variable
        self.subjets    = self.build_dict(cfg.ana_vars['sjvars'])
        self.genfatjets = self.build_dict(cfg.ana_vars['genak8jets'])
        self.gensubjets = self.build_dict(cfg.ana_vars['gensubjets'])
        ak8jmsjmr_helper(self, self.jec_sys) # apply corrections to softdrop mass
        self.fatjets = self.tmp_fatjets  # hand over correct ak8 info
        if __name__ == "__main__": # if doing testing 
            self.fatjets['FatJet_msoftdrop'] = self.tmp_fatjets["FatJet_msoftdrop_altnosmear"]
        del self.subjets, self.genfatjets, self.gensubjets, self.tmp_fatjets
        ## end of in-house JMR
        #
        self.jets      = self.build_dict(
            cfg.ana_vars['ak4vars']+cfg.ana_vars['ak4lvec']['TLVars']+(
                [] if self.isData else (cfg.ana_vars['ak4mcvars'] if self.jec_sys is None else ['Jet_btagSF_deepcsv_shape'])
            ))
        self.electrons = self.build_dict(cfg.lep_sel_vars['electron']) 
        self.muons     = self.build_dict(cfg.lep_sel_vars['muon']) 
        self.events    = self.build_dict(cfg.ana_vars['event']+(
            (cfg.ana_vars['sysvars_mc']+cfg.ana_vars[f'sysvars_{self.year}']) 
            if not self.isData else []))  
        #print(len(self.events['genWeight']))
        # other things like gen info
        if not self.isData:
            self.geninfo    = self.build_dict(cfg.ana_vars['genpvars']) 
            self.lheweights = self.build_dict(cfg.ana_vars['lheWeights'])
            
        self.hlt        = self.build_dict(cfg.ana_vars['dataHLT_all']+cfg.ana_vars[f'dataHLT_{self.year}'])
        #self.subjets    = self.build_dict(cfg.ana_vars['ak8sj'])
        # wont keep
        self.filters    = self.build_dict(cfg.ana_vars['filters_all']+cfg.ana_vars['filters_year'][self.year]) 
        #del self.precut # need later i think
        self.f.close()
        # ===================== #
        # apply object criteria
        self.soft_electrons = self.electrons[self.is_a_soft_electron()]
        self.electrons      = self.electrons[self.is_a_electron()]
        self.soft_muons     = self.muons[    self.is_a_soft_muon()]
        self.muons          = self.muons[    self.is_a_muon()]
        self.jets           = self.jets[     self.is_a_jet()] # still need to lepton clean
        self.fatjets        = self.fatjets[  self.is_a_fatjet()]
        # define event selection
        self.event_mask   = self.get_event_selection()
        # apply event selection
        self.jets            = self.jets[          self.event_mask]
        self.fatjets         = self.fatjets[       self.event_mask]
        self.electrons       = self.electrons[     self.event_mask]
        self.soft_electrons  = self.soft_electrons[self.event_mask]
        self.muons           = self.muons[         self.event_mask]
        self.soft_muons      = self.soft_muons[    self.event_mask]
        self.events          = self.events[        self.event_mask]
        if not self.isData:
            self.geninfo    = self.geninfo[   self.event_mask]
            self.lheweights = self.lheweights[self.event_mask]
        self.hlt        = self.hlt[       self.event_mask]
        # might drop subjets
        #self.subjets    = self.subjets[   self.event_mask]
        #print(len(self.events['genWeight']))
        ''' 
        add interesting info to events: 
        (lepton pt, eta, phi , etc...) 
        njets, nfjets, nbjets
        ps, sc, eft
        hem veto weight
        '''
        self.handle_multiplicity_HEM_info()
        if not self.isData:
            self.handle_lheweights()
            self.Meta.add_btagweightsf_counts(self.jets, self.events, self.geninfo)
        self.handle_lep_info()
        #
        self.btag_event_mask = self.get_b_tag_eventmask()
        self.jets            = self.jets[          self.btag_event_mask]
        self.fatjets         = self.fatjets[       self.btag_event_mask]
        #self.subjets         = self.subjets[       self.btag_event_mask]
        self.soft_electrons  = self.soft_electrons[self.btag_event_mask]
        self.soft_muons      = self.soft_muons[    self.btag_event_mask]
        self.hlt             = self.hlt[           self.btag_event_mask]
        self.events          = self.events [       self.btag_event_mask]
        if not self.isData:
            self.geninfo = self.geninfo[  self.btag_event_mask]
            # add full list of weights to metaData and apply cuts
            pdf_helper.apply_cuts(self)
            pdf_helper.add_pdfweights_to_events(self)
        #
    #
    def get_skim(self):
        __out_dict = {
            'ak4':self.jets,
            'ak8':self.fatjets,
            #'ak8':{**self.fatjets,**self.subjets},
            'soft_e':self.soft_electrons,
            'soft_mu':self.soft_muons,
            #'sj': self.subjets
        }
        events = {**self.events,**self.hlt}
        __out_dict['events'] = pd.DataFrame.from_dict({k:v for k,v in events.items()})
        if not self.isData:
            __out_dict['gen'] = self.geninfo
            __out_dict['metaData'] = self.Meta.get_metadata()
        # close root file
        #self.f.close()
        #
        return __out_dict
            

    # === functions to add info to events === #
    def handle_lheweights(self):
        ps_w  = self.lheweights['PSWeight'].pad(4).fillna(1)
        try:
            sc_w  = self.lheweights['LHEScaleWeight'].pad(9).fillna(1)
        except:
            sc_w = np.ones(shape=(len(ps_w),9))

        self.events['ISR_Up']   = ps_w[:,2]
        self.events['ISR_Down'] = ps_w[:,0]
        self.events['FSR_Up']   = ps_w[:,3]
        self.events['FSR_Down'] = ps_w[:,1]
        #
        self.events['mu_r_Up']    = sc_w[:,7]
        self.events['mu_r_Down']  = sc_w[:,1]
        self.events['mu_f_Up']    = sc_w[:,5]
        self.events['mu_f_Down']  = sc_w[:,3]
        self.events['mu_rf_Up']   = sc_w[:,8]
        self.events['mu_rf_Down'] = sc_w[:,0]
        #
        if 'EFT' in self.sample:
            eft_w =  self.lheweights['LHEReweightingWeight'].pad(184).fillna(1)
            for i in range(184):
                self.events[f'EFT{i}'] = eft_w[:,i]

        
    #
    def handle_lep_info(self): # assumes event filter already applied (1 lepton per event)
        get_lep_info = (lambda k : aj.concatenate(
            [self.muons[f'Muon_{k}'], self.electrons[f'Electron_{k}']], axis=1).flatten()
        )
        single_mu = (self.muons['Muon_pt'].counts         == 1)
        single_el = (self.electrons['Electron_pt'].counts == 1)
        #print(len(get_lep_info('pt')))
        self.events.update({
            'Lep_pt'  : get_lep_info('pt'),
            'Lep_eta' : get_lep_info('eta'),
            'Lep_phi' : get_lep_info('phi'),
            'Lep_mass': get_lep_info('mass'),
            'Lep_iso' : get_lep_info('miniPFRelIso_all'),
            'Lep_ch'  : get_lep_info('charge'),
            'passSingleLepMu'   : single_mu,
            'passSingleLepElec' : single_el
        })
        self.events['Muon_mediumPromptId'] = self.muons['Muon_mediumPromptId'].pad(1).fillna(0)
        #self.events['Muon_tightId']        = self.muons['Muon_tightId'].pad(1).fillna(0)
        self.events['Muon_sip3d']          = self.muons['Muon_sip3d'].pad(1).fillna(999)
        #self.events['Muon_dxy']            = abs(self.muons['Muon_dxy'].pad(1).fillna(999))
        #self.events['Muon_dz']             = abs(self.muons['Muon_dz'].pad(1).fillna(999))
        #
        self.events['Electron_sip3d']         = self.electrons['Electron_sip3d'].pad(1).fillna(999)
    #
    def handle_multiplicity_HEM_info(self):
        self.events.update({
            'n_ak4jets': self.jets['Jet_pt'].counts,
            'n_ak8jets': self.fatjets['FatJet_pt'].counts,
            'nBottoms' : self.jets['Jet_pt'][(self.jets['Jet_btagDeepB'] > cfg.ZHbb_btagWP[self.year])].counts
        })
        if self.year == '2018':
            if not self.isData: # is MC 
                self.events.update({
                    'HEM_weight': np.where(self.get_HEM_veto(), 1, 0.3518) # fraction of lumi with no HEM issue
                })
            else : # isData
                self.events.update({
                    'PassHEM_veto' : np.where(self.events['run']> 319077, self.get_HEM_veto(), True)
                })
        

    # === object criteria functions === #
    @staticmethod
    def is_lep_cleaned(lep, l_k, jets, j_k, cut=0.4): # cut should be 0.4, 0.8 
        # cleans any jets matched with analysis lepton
        lep_eta, lep_phi = lep[f'{l_k}_eta'], lep[f'{l_k}_phi']
        jets_eta, jets_phi = jets[f'{j_k}_eta'], jets[f'{j_k}_phi']
        jet_mask = ak_crosscleaned(lep_eta,lep_phi,jets_eta,jets_phi,cut)
        return jet_mask
        

    def is_a_jet(self):
        return  (
            (self.jets['Jet_pt']       > 30) & 
            (abs(self.jets['Jet_eta']) < 2.4) & 
            ((self.jets['Jet_pt'] > 50) | (self.jets['Jet_puId'] >= 4) ) &
            ( self.jets['Jet_jetId'] >= 2) & 
            ( self.is_lep_cleaned(self.electrons,'Electron',self.jets,'Jet',0.4) == True ) &  
            ( self.is_lep_cleaned(self.muons,'Muon',self.jets,'Jet',0.4) == True ) 
        )
    #
    def is_a_fatjet(self):
        return (
            (self.fatjets['FatJet_pt'] >  200) &
            (abs(self.fatjets['FatJet_eta']) < 2.4) &   
            (self.fatjets['FatJet_msoftdrop'] >= 50) & 
            (self.fatjets['FatJet_msoftdrop'] <= 200) & 
            ( self.fatjets['FatJet_jetId'] >= 2) &
            ( self.is_lep_cleaned(self.electrons,'Electron',self.fatjets,'FatJet',0.8) == True ) & 
            ( self.is_lep_cleaned(self.muons,'Muon',self.fatjets,'FatJet',0.8) == True ) 
        )
                        
    #
    def is_a_electron(self):
        return  cfg.lep_sel['electron'][self.year](self.electrons)
    def is_a_soft_electron(self):
        return  (
            (self.electrons['Electron_cutBasedNoIso'] >= 1) &
            (abs(self.electrons['Electron_eta']) < 2.5) &
            (self.electrons['Electron_pt'] <= (30 if self.year == '2016' else 35))
        )
    #
    def is_a_muon(self):
        return cfg.lep_sel['muon'](self.muons)
    def is_a_soft_muon(self):
        return  (
            ((self.muons['Muon_softId'] == 1) | (self.muons['Muon_mediumId'] == 1)) &
            (self.muons['Muon_pt'] <= 30) &
            (abs(self.muons['Muon_eta']) < 2.4)
        )
    # === event criteria functions === #
    def get_MET_filter(self) :    
        return ((self.filters['Flag_goodVertices']                       == 1)       &
                (self.filters['Flag_globalSuperTightHalo2016Filter']     == 1)       &
                (self.filters['Flag_HBHENoiseFilter']                    == 1)       & 
                (self.filters['Flag_HBHENoiseIsoFilter']                 == 1)       & 
                (self.filters['Flag_EcalDeadCellTriggerPrimitiveFilter'] == 1)       & 
                (self.filters['Flag_BadPFMuonFilter']                    == 1)       &
                ((self.filters['Flag_eeBadScFilter'] == 1) if self.isData else True) &
                ((self.filters['Flag_ecalBadCalibFilterV2'] == 1) 
                 if (self.year == '2017' or self.year == '2018') else True )
        )
    #
    def get_HEM_veto(self) : 
        elec_hem = lambda : ((self.electrons['Electron_pt'] > 20 )    &
                             (self.electrons['Electron_eta'] > -3.0)  &
                             (self.electrons['Electron_eta'] < -1.4)  &
                             (self.electrons['Electron_phi'] < -0.87) &
                             (self.electrons['Electron_phi'] > -1.57))
        #returns True if not in problematic region
        return  (self.electrons['Electron_pt'][elec_hem()].counts > 0) == False 
        
    #
    def pass_goldenjson(self):
        run , lumi = self.events['run'].flatten(), self.events['luminosityBlock'].flatten()
        pass_runlumi = np.zeros_like(run)
        unique_rl_pairs = np.unique(np.column_stack((run,lumi)),axis=0)
        for r,l in unique_rl_pairs:
            if str(r) in self.golden_json:
                for start,finish in self.golden_json[str(r)]:
                    if l in range(start,finish+1): # add one to make range effectively inclusive
                        #print(r,self.events['run'] == r)
                        #print(l,self.events['luminosityBlock'] == l)
                        pass_runlumi = np.where(
                            ((self.events['run'] == r) & (self.events['luminosityBlock'] == l)).flatten(),1,pass_runlumi)
        return pass_runlumi
        

    def get_event_selection(self): # after objects have been defined
        return ( (self.jets['Jet_pt'].counts >= (5 if self.jec_sys is None else 5)) & # n_jets >= 5 is the norm
                 (self.fatjets['FatJet_pt'].counts >= 1) &
                 (self.events['MET_pt'] > 20) &
                 (self.electrons['Electron_pt'].counts + self.muons['Muon_pt'].counts == 1) &
                 self.get_MET_filter() &
                 (self.pass_goldenjson() == 1 if self.isData else True) &
                 (
                     np.where(self.events['run'] >= 319077, self.get_HEM_veto(), True) if 
                     (self.year == '2018' and self.isData) else True
                 )
             )
    def get_b_tag_eventmask(self):
        return (self.events['nBottoms'] >= 2)
    # === helper functions === #
    def set_tree_from_roofile(self,roofile, estart=None, estop=None):
        ''' Set tree array method and tree pandas method '''
        self.f = uproot.open(self.roofile)
        #print(estart, estop)
        tree         =  self.f.get('Events') # hardcoded
        self.tarray  = self.tarray_wrapper(functools.partial(tree.array,      entrystart=estart, entrystop=estop))
        self.tpandas = functools.partial(tree.pandas.df , entrystart=estart, entrystop=estop)
        self.precut  = self.set_pre_cuts()
        return tree

    def set_pre_cuts(self):
        _c_keys = ['nJet','nFatJet','MET_pt','nMuon','nElectron']
        executor = concurrent.futures.ThreadPoolExecutor()
        _c_vars = AnaDict({k: self.tarray(k, executor=executor) for k in _c_keys})
        _precut = ( (_c_vars['nMuon'] + _c_vars['nElectron'] >= 1) &
                    (_c_vars['nJet'] >= 5)  &
                    (_c_vars['nFatJet'] >= 1) &
                    (_c_vars['MET_pt'] >= 20) )
        return _precut
        

    @staticmethod
    def tarray_wrapper(tarray):
        def wrapper(key, **kwargs):
            try:
                out = tarray(key, **kwargs)
            except KeyError:
                print(f"{key} is not in root file!!! Returning zeros")
                out = np.zeros_like(tarray('run',**kwargs).flatten())
            return out
        return wrapper
        
    def build_dict(self, keys, with_interp=True):
        executor = concurrent.futures.ThreadPoolExecutor()
        if with_interp:
            return AnaDict({k: self.tarray(self.ana_vars[k], executor=executor)[self.precut] for k in keys})
            #return AnaDict({k: self.tarray(self.ana_vars[k], executor=executor) for k in keys})
        else:
            return AnaDict({k: self.tarray(k, executor=executor)[self.precut] for k in keys})
            #return AnaDict({k: self.tarray(k, executor=executor) for k in keys})
    # === ~ Skim class === #
    # extra local test stuff
    def local_test(self):
        # write to file to compare with Andrew' Skimmer
        exit()
        out_txt_file = open(f"bryan_{self.year}.txt", "w")
        for i,j,k in zip(self.events['luminosityBlock'],self.events['event'],self.events['run']):
            print(f"Run {k}, LS {i}, event {j}")
            out_txt_file.writelines([f"Run {k}, LS {i}, event {j}\n"])
        out_txt_file.close()
        #print(self.Meta.get_metadata())
        print(f"{len(self.events['event'])} events passed")
        print(self.events[self.events['event']==4644390569])
        #print(self.events[ self.events['event']==804875471]['Lep_eta'],
        #      self.events[ self.events['event']==804875471]['Lep_phi'])
                
if __name__ == '__main__':
    #test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2018/ttHTobb_2018/839BA380-7826-9140-8C16-C5C0903EE949_Skim_12.root'
    #sample='TTToSemiLeptonic'
    #test_file  = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2018/TTToSemiLeptonic_2018/D6501B6C-8B76-BF42-B677-64680733A780_Skim_19.root'
    #test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2017/TTToSemiLeptonic_2017/B8B652A8-ED88-2F48-9979-637E36F30138_Skim_72.root'
    #test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2017/TTToSemiLeptonic_2017/DEDD55D3-8B36-3342-8531-0F2F4C462084_Skim_134.root' 
    #test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2016/TTToSemiLeptonic_2016/CA4521C3-F903-8E44-93A8-28F5D3B8C5E8_Skim_121.root'
    #test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2016/ttHTobb_2016/A1490EBE-FA8A-DE40-97F8-FCFBAB716512_Skim_11.root'
    #sample = "ttHTobb"
    #sample = 'TTbb_2L2Nu'
    #test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed//2017/TTbb_2L2Nu_2017/prod2017MC_v7_NANO_2_Skim_11.root'
    #sample= 'TTZToBB'
    #test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2016/TTZToBB_2016/CA01B0AA-229F-E446-B4FE-9F2EA2969FAB_Skim_2.root'
    sample = 'ttHToNonbb'
    test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed//2016/ttHToNonbb_2016/6D0E1ED8-6F95-FE45-8967-96805FBF1818_Skim_25.root'
    #sample = 'Data_SingleMuon'
    #test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2017/Data_SingleMuon_2017_PeriodD/9A53E75E-4E0E-5A4A-A8C3-A91333DA906D_Skim_8.root'
    #test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2018/Data_SingleMuon_2018_PeriodC/C8A1B18B-F06D-7D4B-80AC-3FD4774625AF_Skim_14.root'
    #test_file = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2016/Data_SingleMuon_2016_PeriodE/E1D9A669-F314-6944-9946-93C893F180A9_Skim_10.root'
    #test_file    = '/eos/uscms/store/user/bcaraway/SingleE_test_2017.root'
    year = re.search(r"201\d", test_file).group()
    golden_json=json.load(open(cfg.goodLumis_file[year]))
    _ = Skim(test_file, sample, year, isData='Data' in sample, jec_sys=None, golden_json=golden_json)
    _.local_test()
    #AnaDict(_.get_skim()).to_pickle('SingleE_2017.pkl')
    
