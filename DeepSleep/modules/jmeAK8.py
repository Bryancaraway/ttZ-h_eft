########################
### Get data         ###
### for analysis     ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#
import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
    print(sys.path[1])
#
import uproot
from coffea.analysis_objects import JaggedCandidateArray, JaggedTLorentzVectorArray
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty, JetTransformer, JetResolution, JetResolutionScaleFactor
from coffea.lookup_tools import extractor
import awkward as ak
#
import config.jme_cff as cfg
from lib.fun_library import fillne, t2Run
from modules.AnaDict import AnaDict
#from modules.AnaVars import AnaVars
#
import numpy as np
np.random.seed(0)
import pandas as pd
##

class JMEAK8 :
    '''
    user coffea module to calculate ak8 related uncertainties
    JES, and JER vairations on pt and sd_Mass
    '''
    year = None
    
    def __init__(self, year):
        self.year = year
        self.jec_info   = cfg.ak8_jme_files['jec'](self.year)
        self.junc_info  = cfg.ak8_jme_files['junc'](self.year)
        self.jer_info   = cfg.ak8_jme_files['jer'](self.year)
        self.jersf_info = cfg.ak8_jme_files['jersf'](self.year)
        #
        fj_ext = extractor()
        fj_ext.add_weight_sets(
            [f"* * {jec_i}" for jec_i in self.jec_info]
            +[f"* * {self.junc_info}",
              f"* * {self.jer_info}",
              f"* * {self.jersf_info}"]
        )
        fj_ext.finalize()
        self.fj_eval = fj_ext.make_evaluator()

    def prepare_transformer(self):
        ''' 
        use JetCorrectionUncertainty and 
        JetResolutionScaleFactor to get ak8 uncertainties 
        '''
        get_name = (lambda x : x.split('/')[-1].split('.j')[0]) # remove directory and .j*.txt from name
        jec_names  = [get_name(jec_i) for jec_i in self.jec_info] # 3 levels = 3 files 
        junc_name  = get_name(self.junc_info) 
        jer_name   = get_name(self.jer_info)
        jersf_name = get_name(self.jersf_info)

        self.jec   = FactorizedJetCorrector(**{n:self.fj_eval[n] for n in jec_names})
        self.junc  = JetCorrectionUncertainty(**{junc_name: self.fj_eval[junc_name]})
        self.jer   = JetResolution(**{jer_name: self.fj_eval[jer_name]})
        self.jersf = JetResolutionScaleFactor(**{jersf_name: self.fj_eval[jersf_name]}) 
        #
        self.Jet_transformer = JetTransformer(jec=self.jec, junc=self.junc, jer=self.jer, jersf=self.jersf)

    def transform_AK8(self, ak8_vars):
        '''
        raw jagged format, from AnaDict
        need to gen match GenAk8Jet to AK8Jet
        '''
        fj = JaggedCandidateArray.candidatesfromcounts(ak8_vars['FatJet_pt_drLeptonCleaned'].counts,
                                                       pt=ak8_vars['FatJet_pt_drLeptonCleaned'].flatten(),
                                                       eta=ak8_vars['FatJet_eta_drLeptonCleaned'].flatten(),
                                                       phi=ak8_vars['FatJet_phi_drLeptonCleaned'].flatten(),
                                                       mass=ak8_vars['FatJet_msoftdrop_drLeptonCleaned'].flatten())

        # have to pad because shitty argmatch function doesnt really like empty arrays
        genfj = JaggedCandidateArray.candidatesfromcounts(ak8_vars['GenJetAK8_pt'].pad(1).fillna(0.1).counts,
                                                          pt=ak8_vars  ['GenJetAK8_pt'].pad(1).fillna(0.1).flatten(),
                                                          eta=ak8_vars ['GenJetAK8_eta'].pad(1).fillna(-5).flatten(),
                                                          phi=ak8_vars ['GenJetAK8_phi'].pad(1).fillna(0).flatten(),
                                                          mass=ak8_vars['GenJetAK8_mass'].pad(1).fillna(0.1).flatten())
        
        mb= fj.match(genfj, deltaRCut=0.4)
        fj.add_attributes(ptGenJet=genfj.pt[fj.argmatch(genfj, deltaRCut=0.4)]) # will retrun -1 index value if no match so have to take care of this
        fj.ptGenJet[~mb] = 0. # have to do this manually for jer to work
        fj.add_attributes(ptRaw = fj.pt * (1 - ak8_vars['FatJet_rawFactor_drLeptonCleaned']),
                          massRaw = fj.mass * (1 - ak8_vars['FatJet_rawFactor_drLeptonCleaned']))
        # thanks jetmet group for wasting my time...
        fj['rho']  = fj.pt.ones_like()   # this actually doesnt really do anything...
        #fj['area'] = fj.pt.ones_like()*2 # this actually doesnt really do anything...
        fj['area'] = ak8_vars['FatJet_area_drLeptonCleaned']
        #
        #old = fj.pt # remove !!!
        self.Jet_transformer.transform(fj)
        # prepare jet collection to return with the right format for my anaylsis
        #'pt_jer_up', 'mass_jer_up', 'pt_jer_down', 'mass_jer_down', 'pt_jes_up', 'mass_jes_up', 'pt_jes_down', 'mass_jes_down']
        _fj = {}
        _interp = {'pt':'pt',
                   'mass':'msoftdrop'}
        for k in ['pt','mass']:
            _fj[f'FatJet_{_interp.get(k)}_drLeptonCleaned']              = fj[f'__fast_{k}']
            _fj[f'FatJet_{_interp.get(k)}_jesTotalUp_drLeptonCleaned']   = fj[f'{k}_jes_up']
            _fj[f'FatJet_{_interp.get(k)}_jesTotalDown_drLeptonCleaned'] = fj[f'{k}_jes_down']
            _fj[f'FatJet_{_interp.get(k)}_jerUp_drLeptonCleaned']        = fj[f'{k}_jer_up']
            _fj[f'FatJet_{_interp.get(k)}_jerDown_drLeptonCleaned']      = fj[f'{k}_jer_down']
        #
        #import matplotlib.pyplot as plt      # remove!!
        #print(abs((old-fj.pt)/old).flatten()[:20])
        #plt.hist(abs((old-fj.mass)/old).flatten(), # remove!!
        #         range=(0,0.1))
        #plt.show()                           # remove!!
        #print(fj.mass)                       # remove!!
        return _fj


if __name__ == '__main__':
    ''' For module testing '''
    x = JMEAK8('2017')
    x.prepare_transformer()
    fj = AnaDict.read_pickle(sys.path[1]+'files/2017/mc_files/TTZH_ak8.pkl')
    #
    fj_array = JaggedCandidateArray.candidatesfromcounts(fj['FatJet_pt_drLeptonCleaned'].counts,
                                                         pt=fj['FatJet_pt_drLeptonCleaned'].flatten(),
                                                         eta=fj['FatJet_eta_drLeptonCleaned'].flatten(),
                                                         phi=fj['FatJet_phi_drLeptonCleaned'].flatten(),
                                                         mass=fj['FatJet_mass_drLeptonCleaned'].flatten())
