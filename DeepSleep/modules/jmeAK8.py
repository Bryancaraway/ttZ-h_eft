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
            [f"* * {j}" for j in self.jec_info]
            +[f"* * {j}" for j in self.junc_info]
            +[f"* * {j}" for j in self.jer_info]
            +[f"* * {j}" for j in self.jersf_info]
        )
        fj_ext.finalize()
        self.fj_eval = fj_ext.make_evaluator()

    def prepare_transformer(self):
        ''' 
        use JetCorrectionUncertainty and 
        JetResolutionScaleFactor to get ak8 uncertainties 
        '''
        get_name = (lambda x : x.split('/')[-1].split('.j')[0]) # remove directory and .j*.txt from name
        jec_names  = [get_name(j) for j in self.jec_info] # 3 levels = 3 files 
        junc_names  = [get_name(j) for j in self.junc_info]
        jer_names   = [get_name(j) for j in self.jer_info]
        jersf_names = [get_name(j) for j in self.jersf_info]

        self.jec_ak8   = FactorizedJetCorrector(  **{n:self.fj_eval[n] for n in jec_names if 'AK4' not in n})
        self.junc_ak8  = JetCorrectionUncertainty(**{n:self.fj_eval[n] for n in junc_names if 'AK4' not in n})
        self.jer_ak8   = JetResolution(           **{n:self.fj_eval[n] for n in jer_names if 'AK4' not in n})
        self.jersf_ak8 = JetResolutionScaleFactor(**{n:self.fj_eval[n] for n in jersf_names if 'AK4' not in n})
        #
        self.Jt_ak8pfpuppi_full = JetTransformer(jec=self.jec_ak8, junc=self.junc_ak8, jer=self.jer_ak8, jersf=self.jersf_ak8)
        self.Jt_ak8pfpuppi_jes  = JetTransformer(jec=self.jec_ak8)

        self.jec_sdm   = FactorizedJetCorrector(  **{n:self.fj_eval[n] for n in jec_names if 'AK8' not in n})
        self.junc_sdm  = JetCorrectionUncertainty(**{n:self.fj_eval[n] for n in junc_names if 'AK8' not in n})
        self.jer_sdm   = JetResolution(           **{n:self.fj_eval[n] for n in jer_names if 'AK8' not in n})
        self.jersf_sdm = JetResolutionScaleFactor(**{n:self.fj_eval[n] for n in jersf_names if 'AK8' not in n})
        #
        self.Jt_ak4pfpuppi_full = JetTransformer(jec=self.jec_sdm, junc=self.junc_sdm, jer=self.jer_sdm, jersf=self.jersf_sdm)
        self.Jt_ak4pfpuppi_jes  = JetTransformer(jec=self.jec_sdm)

    def transform_AK8(self, ak8_vars):
        '''
        raw jagged format, from AnaDict
        need to gen match GenAk8Jet to AK8Jet
        '''
        self.fj = JaggedCandidateArray.candidatesfromcounts(ak8_vars['FatJet_pt_drLeptonCleaned'].counts,
                                                       pt=ak8_vars['FatJet_pt_drLeptonCleaned'].flatten(),
                                                       eta=ak8_vars['FatJet_eta_drLeptonCleaned'].flatten(),
                                                       phi=ak8_vars['FatJet_phi_drLeptonCleaned'].flatten(),
                                                       mass=ak8_vars['FatJet_msoftdrop_drLeptonCleaned'].flatten())
        print(self.fj.pt)
        # have to pad because shitty argmatch function doesnt really like empty arrays
        genfj = JaggedCandidateArray.candidatesfromcounts(ak8_vars['GenJetAK8_pt'].pad(1).fillna(0.1).counts,
                                                          pt=ak8_vars  ['GenJetAK8_pt'].pad(1).fillna(0.1).flatten(),
                                                          eta=ak8_vars ['GenJetAK8_eta'].pad(1).fillna(-5).flatten(),
                                                          phi=ak8_vars ['GenJetAK8_phi'].pad(1).fillna(0).flatten(),
                                                          mass=ak8_vars['GenJetAK8_mass'].pad(1).fillna(0.1).flatten())
        
        mb= self.fj.match(genfj, deltaRCut=0.8)
        self.fj.add_attributes(ptGenJet=genfj.pt[self.fj.argmatch(genfj, deltaRCut=0.8)]) # will retrun -1 index value if no match so have to take care of this
        del genfj
        
        self.fj.ptGenJet[~mb] = 0. # have to do this manually for jer to work
        self.fj.add_attributes(ptRaw = self.fj.pt * (1 - ak8_vars['FatJet_rawFactor_drLeptonCleaned']),
                               massRaw = self.fj.mass * (1 - ak8_vars['FatJet_rawFactor_drLeptonCleaned']))
        # thanks jetmet group for wasting my time...
        self.fj['rho']  = self.fj.pt.ones_like()   # this actually doesnt really do anything...
        #fj['area'] = fj.pt.ones_like()*2 # this actually doesnt really do anything...
        self.fj['area'] = ak8_vars['FatJet_area_drLeptonCleaned']
        #
        #old = fj.pt # remove !!!
        self.Jt_ak8pfpuppi_full.transform(self.fj)
        # prepare jet collection to return with the right format for my anaylsis
        #'pt_jer_up', 'mass_jer_up', 'pt_jer_down', 'mass_jer_down', 'pt_jes_up', 'mass_jes_up', 'pt_jes_down', 'mass_jes_down']
        _fj = {}
        print(self.fj.pt)
        print(self.fj[f'__fast_pt'])
        for k in ['pt']:#,'mass']:
            #_fj[f'FatJet_{k}_drLeptonCleaned']              = self.fj[f'__fast_{k}']  
            _fj[f'FatJet_{k}_jesTotalUp_drLeptonCleaned']   = self.fj[f'{k}_jes_up']  
            _fj[f'FatJet_{k}_jesTotalDown_drLeptonCleaned'] = self.fj[f'{k}_jes_down']
            _fj[f'FatJet_{k}_jerUp_drLeptonCleaned']        = self.fj[f'{k}_jer_up']  
            _fj[f'FatJet_{k}_jerDown_drLeptonCleaned']      = self.fj[f'{k}_jer_down']
        #
        _fj['FatJet_msoftdrop_drLeptonCleaned']              = ak8_vars['FatJet_msoftdrop_drLeptonCleaned'] * self.fj['__fast_pt']  / self.fj.pt
        _fj['FatJet_msoftdrop_jesTotalUp_drLeptonCleaned']   = ak8_vars['FatJet_msoftdrop_drLeptonCleaned'] * self.fj['pt_jes_up']  / self.fj.pt
        _fj['FatJet_msoftdrop_jesTotalDown_drLeptonCleaned'] = ak8_vars['FatJet_msoftdrop_drLeptonCleaned'] * self.fj['pt_jes_down']/ self.fj.pt
        _fj['FatJet_msoftdrop_jerUp_drLeptonCleaned']        = ak8_vars['FatJet_msoftdrop_drLeptonCleaned'] * self.fj['pt_jer_up']  / self.fj.pt
        _fj['FatJet_msoftdrop_jerDown_drLeptonCleaned']      = ak8_vars['FatJet_msoftdrop_drLeptonCleaned'] * self.fj['pt_jer_down']/ self.fj.pt

        #import matplotlib.pyplot as plt
        #plt.hist(_fj['FatJet_msoftdrop_drLeptonCleaned'].flatten(), range=(50,200), histtype='step', bins=[50,80,105,145,200], label='nom')
        #plt.hist(_fj['FatJet_msoftdrop_jesTotalUp_drLeptonCleaned'].flatten(), range=(50,200), histtype='step', bins=[50,80,105,145,200], label='jesup')
        #plt.hist(_fj['FatJet_msoftdrop_jesTotalDown_drLeptonCleaned'].flatten(), range=(50,200), histtype='step', bins=[50,80,105,145,200], label='jesdown')
        #plt.hist(_fj['FatJet_msoftdrop_jerUp_drLeptonCleaned'].flatten(), range=(50,200), histtype='step', bins=[50,80,105,145,200], label='jerup')
        #plt.hist(_fj['FatJet_msoftdrop_jerDown_drLeptonCleaned'].flatten(), range=(50,200), histtype='step', bins=[50,80,105,145,200], label='jerdown')
        #plt.legend()
        #plt.show()
        #exit()
        #
        return _fj
                   
    def transform_SDM(self, ak8_vars):
        # assemble 4 vector per ak8 jet made up of 2 sj
        sj1, sj2 = ak8_vars['FatJet_subJetIdx1_drLeptonCleaned'], ak8_vars['FatJet_subJetIdx2_drLeptonCleaned']
        #sj1 = sj1[sj1 != -1] # prepare for indexing 
        #sj2 = sj2[sj2 != -1]
        # assemble subjet array and gensubjet array
        sj = JaggedCandidateArray.candidatesfromcounts(ak8_vars['SubJet_pt'].pad(1).fillna(0.1).counts,
                                                       pt=  (ak8_vars['SubJet_pt']* (1-ak8_vars['SubJet_rawFactor'])).pad(1).fillna(0.1).flatten(),
                                                       eta= ak8_vars['SubJet_eta'].pad(1).fillna(-99).flatten(),
                                                       phi= ak8_vars['SubJet_phi'].pad(1).fillna(0).flatten(),
                                                       mass=(ak8_vars['SubJet_mass']* (1-ak8_vars['SubJet_rawFactor'])).pad(1).fillna(0.1).flatten())

        gensj = JaggedCandidateArray.candidatesfromcounts(ak8_vars['SubGenJetAK8_pt'].pad(1).fillna(0.1).counts,
                                                          pt=ak8_vars  ['SubGenJetAK8_pt'].pad(1).fillna(0.1).flatten(),
                                                          eta=ak8_vars ['SubGenJetAK8_eta'].pad(1).fillna(-5).flatten(),
                                                          phi=ak8_vars ['SubGenJetAK8_phi'].pad(1).fillna(0).flatten(),
                                                          mass=ak8_vars['SubGenJetAK8_mass'].pad(1).fillna(0.1).flatten())
        mb= sj.match(gensj, deltaRCut=0.4)
        argm = sj.argmatch(gensj, deltaRCut=0.4)
        # okay so lets correct the subjets 
        sj.add_attributes(ptGenJet=gensj.pt,
                          ptRaw = sj.pt,
                          massRaw = sj.mass)
        sj['rho'] = sj.pt.ones_like()
        sj['area'] = sj.pt.ones_like()
        #self.Jt_ak4pfpuppi_full.transform(sj) 
        #
        gensj = gensj[argm]
        gensj.pt[~mb] = 0. # have to do this manually for jer to work
        gensj.eta[~mb] = -5.
        gensj.phi[~mb] = 0
        gensj.mass[~mb] = 0.

        #
        groomed = sj[sj1].p4 + sj[sj2].p4
        gengroomed = gensj[sj1].p4 + gensj[sj2].p4
        del gensj, sj
        #
        #self.Jt_ak4pfpuppi_jes.transform(self.fj)
        #fjpt_corr = self.fj.pt # will need to also keep jes up/donw
        #self.Jt_ak8pfpuppi_full.transform(self.fj)
        fjpt_smear = self.fj.pt # will need to also keep jes up/down
        #fjpt_jerNomVal = fjpt_smear/fjpt_corr 

        puppi_corr = cfg.puppicorr_gen[self.year](fjpt_smear.flatten())*np.where(abs(self.fj.eta.flatten()) <= 1.3, 
                                                                                cfg.puppicorr_reco0eta1p3[self.year](fjpt_smear.flatten()), 
                                                                                cfg.puppicorr_reco1p3eta2p5[self.year](fjpt_smear.flatten()))
        puppi_sdm_corr = ak.JaggedArray.fromoffsets(self.fj.pt.offsets, puppi_corr)
        # get jet mass scale
        puppi_sdm_jms  = cfg.puppi_sdm_jms_jmr['jms'][self.year]
        # corrct sdm with puppi corrections
        groomed_corr = JaggedCandidateArray.candidatesfromcounts(groomed.counts,
                                                                 pt=groomed.pt.flatten(),
                                                                 eta=groomed.eta.flatten(), 
                                                                 phi=groomed.phi.flatten(),
                                                                 mass=(groomed.mass * puppi_sdm_corr).flatten())
        #fj_sdm_corr = groomed.mass * puppi_sdm_corr
        # get jet mass smear
        puppi_sdm_jmrNomVal = ak.JaggedArray.fromoffsets(self.fj.pt.offsets, self.get_jmr(groomed.pt,groomed.eta,groomed.mass, gengroomed.mass))
        # assemble softdropmass fully corrected with uncertainties
        #fj_sdm = fjpt_jerNomVal*puppi_sdm_jms*puppi_sdm_jmrNomVal*groomed_corr.mass
        fj_sdm = groomed_corr.mass#*puppi_sdm_jms*puppi_sdm_jmrNomVal
        

        # finally get rid of sdm where it was originally -1 
        fj_sdm = ak.JaggedArray.fromoffsets(self.fj.pt.offsets, np.where(ak8_vars['FatJet_msoftdrop_drLeptonCleaned'].flatten() == -1,
                                                                         -1,
                                                                         fj_sdm.flatten()))
        print(fj_sdm)
        print(ak8_vars['FatJet_msoftdrop_drLeptonCleaned'])
        #
        import matplotlib.pyplot as plt
        plt.hist(fj_sdm.flatten()[self.fj.pt.flatten() > 300], range=(50,200), histtype='step', bins=[50,80,105,145,200], label='new_corr')
        plt.hist(ak8_vars['FatJet_msoftdrop_drLeptonCleaned'].flatten()[self.fj.pt.flatten() > 300], range=(50,200), histtype='step', bins=[50,80,105,145,200], label='old_corr')
        plt.hist(gengroomed.mass.flatten()[self.fj.pt.flatten() > 300], range=(50,200), histtype='step', bins=[50,80,105,145,200], label='gen')
        plt.legend()
        plt.show()
        # need to reconstruct the (gen)groomed jet by added the two (gen)sj together
        exit()
        
    def get_jmr(self,gpt, geta, sdm, gen_sdm):
        puppi_sdm_jmr  = cfg.puppi_sdm_jms_jmr['jmr'][self.year]
        jer = np.where(abs(geta.flatten()) <= 1.3, 
                       cfg.puppicorr_massReso_0eta1p3(gpt.flatten()),
                       cfg.puppicorr_massReso_1p3eta2p5(gpt.flatten()))
        
        jersmear = jer*np.random.normal(size=jer.size)
        doHybrid = gen_sdm.flatten() > 0
        jsmear = np.where(doHybrid,
                          1 + (puppi_sdm_jmr - 1) * (sdm - gen_sdm).flatten() / (sdm.flatten()),
                          1. + np.sqrt(np.maximum(puppi_sdm_jmr**2 - 1.0, 0)) * jersmear
                      )
        return np.where(jsmear*sdm.flatten() < 1.e-2, 1.e-2, jsmear)


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
