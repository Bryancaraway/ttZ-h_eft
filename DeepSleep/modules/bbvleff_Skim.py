########################
### Skim data        ###
### for AK8 bbvL     ###
### efficiecny sf    ###
### calc             ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#
import os
import sys
import json
import numpy as np
import pandas as pd
#from modules.metaSkim import SkimMeta
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import config.ana_cff as cfg
from lib.fun_library import argmatch, t2Run, deltaR, fillne, save_pdf
from modules.AnaDict import AnaDict
from modules.AnaVars import AnaVars
from modules.metaSkim import SkimMeta
from modules.ak8jmsjmr_helper import ak8jmsjmr_helper
from modules.pdfweight_helper import PDFHelper
from modules.Skim import Skim


class bbvLEffSkim(Skim) :
    '''
    This class borrows the functions
    and obj selection definitions from
    the super class: SKim
    --> with changes to event selection
    reqs, and only fatjet reco info
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)

    def startSkim(self):
        self.ana_vars = AnaVars(self.year, self.isData, jec_sys=self.jec_sys) 
        self.tree     = self.set_tree_from_roofile(self.roofile)
        # handle pdfweight calc/matching if needed
        #if not self.isData:
        pdf_helper = PDFHelper(self) # currently does nothing for EFT samples
        # prepMeta metadata factory
        #self.Meta = SkimMeta(self.sample, self.year, self.isData, self.tree, self.jec_sys, pdf=pdf_helper.pdfweights)
        # define event information
        # apply in house jmr 
        self.fatjets   = self.build_dict(
            cfg.ana_vars['ak8lvec']['TLVars']+['FatJet_jetId','FatJet_deepTagMD_bbvsLight','FatJet_msoftdrop'], with_interp=False)
        self.jets      = self.build_dict(cfg.ana_vars['ak4vars']+cfg.ana_vars['ak4lvec']['TLVars'])
        self.events    = self.build_dict(cfg.ana_vars['event']+(
            (cfg.ana_vars['sysvars_mc']+cfg.ana_vars[f'sysvars_{self.year}']) 
            if not self.isData else []))  
        #print(len(self.events['genWeight']))
        # other things like gen info
        if not self.isData:
            self.geninfo    = self.build_dict(['GenPart_eta', 'GenPart_phi', 
                                               'GenPart_pdgId', 'GenPart_genPartIdxMother', 'GenPart_status']) 
            #self.geninfo    = self.geninfo[((abs(self.geninfo['GenPart_pdgId']) <= 6) |
            #                                ((abs(self.geninfo['GenPart_pdgId']) >= 20) & 
            #                                 (abs(self.geninfo['GenPart_pdgId']) <= 25)))]
    
            
        #self.subjets    = self.build_dict(cfg.ana_vars['ak8sj'])
        # wont keep
        self.filters    = self.build_dict(cfg.ana_vars['filters_all']+cfg.ana_vars['filters_year'][self.year]) 
        #del self.precut # need later i think
        self.f.close()
        # ===================== #
        # apply object criteria
        self.fatjets        = self.fatjets[  self.is_a_fatjet()]
        self.jets        = self.jets[  self.is_a_jet()]
        # define event selection
        self.event_mask   = self.get_event_selection()
        # apply event selection
        self.fatjets         = self.fatjets[       self.event_mask]
        self.jets         = self.jets[       self.event_mask]
        self.events          = self.events[        self.event_mask]
        if not self.isData:
            self.geninfo    = self.geninfo[   self.event_mask]
        # now calculate the bbvl efficiency
        self.calc_eff()
    
    def calc_eff(self):
        g, gen = 'GenPart_', self.geninfo
        fj  = self.fatjets
        if  'QCD'  in self.sample: # measure the fake rate
            gencut = ((abs(gen[g+'pdgId']) <= 6) | (gen[g+'pdgId'] == 21))
            del self.jets
        else: # this is for measuring the signal efficiecny
            print(gen[g+'pdgId'])
            isbb_fromZH     = ((abs(gen[g+'pdgId']) == 5) & ((gen[g+'pdgId'][gen[g+'genPartIdxMother']] == 23) | 
                                                         (gen[g+'pdgId'][gen[g+'genPartIdxMother']] == 25)))
            ispart_fromWfromtop = (((abs(gen[g+'pdgId']) == 11) | (abs(gen[g+'pdgId']) == 13) | 
                                    (abs(gen[g+'pdgId']) == 15) | (abs(gen[g+'pdgId']) < 6)) & 
                                   ((abs(gen[g+'pdgId'][gen[g+'genPartIdxMother'][gen[g+'genPartIdxMother']]]) == 6) & 
                                    (abs(gen[g+'pdgId'][gen[g+'genPartIdxMother']]) ==24)))
            #
            isbb_fromtt     = ((abs(gen[g+'pdgId']) == 5) & (abs(gen[g+'pdgId'][gen[g+'genPartIdxMother']]) == 6))
            ispart_fromtt = (isbb_fromtt | ispart_fromWfromtop)
            ispart_fromisr = ( ((gen[g+'pdgId'] == 21) | (abs(gen[g+'pdgId']) < 6)) & (gen[g+'status'] >= 41) & (gen[g+'status'] <= 49))
            ispart_fromfsr = ( ((gen[g+'pdgId'] == 21) | (abs(gen[g+'pdgId']) < 6)) & (gen[g+'status'] >= 51) & (gen[g+'status'] <= 59))
            #ispart_fromelse = ((ispart_fromtt==False) & (isbb_fromZH==False) & (gen[g+'pdgId'] <= 21) & ((gen[g+'pdgId'] == 21) | (abs(gen[g+'pdgId']) < 6)) & (gen[g+'status'] != 71) )
            #
            gencut = (((gen[g+'pdgId'] == 23) | (gen[g+'pdgId'] == 25)) & (isbb_fromZH.sum() == 2))
        
        fj_gjIdx   = argmatch(fj['FatJet_eta'], fj['FatJet_phi'],
                              gen[g+'eta'][gencut], gen[g+'phi'][gencut],
                              0.6) >= 0
        #if 'QCD' in self.sample:
        genmatched = (fj_gjIdx.sum()==1)
        #matched_fj = fj[(fj_gjIdx) & (fj_gjIdx.sum()==1)].flatten() # numerator and denominator
        #else:
        #    fj_genbbIdx1   = argmatch(fj['FatJet_eta'], fj['FatJet_phi'],
        #                              gen[g+'eta'][(isbb_fromZH & (isbb_fromZH.sum() == 2))], gen[g+'phi'][(isbb_fromZH & (isbb_fromZH.sum() == 2))],
        #                              0.6) >= 0 
        #    fj_genbbIdx2   = argmatch(fj['FatJet_eta'], fj['FatJet_phi'],
        #                              gen[g+'eta'][(isbb_fromZH & (isbb_fromZH.sum() == 2))], gen[g+'phi'][(isbb_fromZH & (isbb_fromZH.sum() == 2))],
        #                              0.6, m_idx=2) >= 0 
        #    both_gbb = (fj_genbbIdx1 & fj_genbbIdx2)
        #    genmatched = ((fj_gjIdx.sum()==1) & (both_gbb.sum()==1)) 
        #    fj_gjIdx = (fj_gjIdx & both_gbb)
        matched_fj = fj[(fj_gjIdx) & (genmatched) ].flatten()
        #matched_fj = fj[(fj_gjIdx) & (fj_gjIdx.sum()==1)] # numerator and denominator
        if 'QCD'  not in self.sample:
            self.jets = self.jets[(genmatched)]
            bb_outfj = np.nansum(deltaR(matched_fj['FatJet_eta'],matched_fj['FatJet_eta'],
                                        fillne(self.jets['Jet_eta'][(self.jets['Jet_btagDeepB'] > cfg.ZHbb_btagWP[self.year])]),
                                        fillne(self.jets['Jet_phi'][(self.jets['Jet_btagDeepB'] > cfg.ZHbb_btagWP[self.year])]),
                                    ) > 0.8, axis=1)
            matched_fj = matched_fj[bb_outfj>=2]
            ispart_fromtt = ispart_fromtt[(genmatched)][bb_outfj>=2]
            gentt           = gen[          (genmatched)][bb_outfj>=2][ispart_fromtt]
            ispart_fromisr = ispart_fromisr[(genmatched)][bb_outfj>=2]
            genisr           = gen[          (genmatched)][bb_outfj>=2][ispart_fromisr]
            ispart_fromfsr = ispart_fromfsr[(genmatched)][bb_outfj>=2]
            genfsr           = gen[          (genmatched)][bb_outfj>=2][ispart_fromfsr]
            #ispart_fromelse = ispart_fromelse[(genmatched)][bb_outfj>=2]
            #genelse         = gen[          (genmatched)][bb_outfj>=2][ispart_fromelse]
            ttpart_matchrZh_1p2 = (deltaR(matched_fj['FatJet_eta'],matched_fj['FatJet_eta'],
                                          fillne(gentt[g+'eta']),fillne(gentt[g+'phi'])) < 1.2)
            ttpart_matchrZh_0p6 = (deltaR(matched_fj['FatJet_eta'],matched_fj['FatJet_eta'],
                                          fillne(gentt[g+'eta']),fillne(gentt[g+'phi'])) < 0.6)
            isrpart_matchrZh_0p6 = (deltaR(matched_fj['FatJet_eta'],matched_fj['FatJet_eta'],
                                          fillne(genisr[g+'eta']),fillne(genisr[g+'phi'])) < 0.6)
            fsrpart_matchrZh_0p6 = (deltaR(matched_fj['FatJet_eta'],matched_fj['FatJet_eta'],
                                          fillne(genfsr[g+'eta']),fillne(genfsr[g+'phi'])) < 0.6)
            print(len(matched_fj['FatJet_pt'][((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) )]))
            print(len(matched_fj['FatJet_pt'][((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) & (np.nansum(isrpart_matchrZh_0p6, axis=1)>0) )]))
            print(len(matched_fj['FatJet_pt'][((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) & (np.nansum(fsrpart_matchrZh_0p6, axis=1)>0) )]))
            print(genfsr[g+'pdgId'][((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) & (np.nansum(fsrpart_matchrZh_0p6, axis=1)>0) )])
            print(genfsr[g+'genPartIdxMother'][((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) & (np.nansum(fsrpart_matchrZh_0p6, axis=1)>0) )])
            print(gen[g+'pdgId'][          (genmatched)][bb_outfj>=2]
                  [((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) & (np.nansum(fsrpart_matchrZh_0p6, axis=1)>0) )][1,:])
            print(gen[g+'genPartIdxMother'][          (genmatched)][bb_outfj>=2]
                  [((matched_fj['FatJet_msoftdrop']>200) & (matched_fj['FatJet_pt']>450) & (np.nansum(ttpart_matchrZh_1p2, axis=1)==0) & (np.nansum(fsrpart_matchrZh_0p6, axis=1)>0) )][1,:])
            exit()
            self.ttpart_matchrZh_0p6 = ttpart_matchrZh_0p6
            self.ttpart_matchrZh_1p2 = ttpart_matchrZh_1p2
            #matched_fj = matched_fj[(np.nansum(ttpart_matchrZh_1p2, axis=1)==0)]
            #ttpart_matchrZh_1p2 = ttpart_matchrZh_1p2[(np.nansum(ttpart_matchrZh_1p2, axis=1)==0)]
        #########
        #ttpart_matchrZh = ttpart_matchrZh[(matched_fj['FatJet_msoftdrop'] > 200) & (matched_fj['FatJet_pt']>700)]
        #partelse_matchrZh = (deltaR(matched_fj['FatJet_eta'],matched_fj['FatJet_eta'],
        #                            fillne(genelse[g+'eta']),fillne(genelse[g+'phi'])) < 1.2)
        #partelse_matchrZh = partelse_matchrZh[(matched_fj['FatJet_msoftdrop'] > 200) & (matched_fj['FatJet_pt']>700)]
        #print(sum(np.nansum(ttpart_matchrZh, axis=1)>0)/len(ttpart_matchrZh), "From matched to top system product")
        #print(sum(np.sum(partelse_matchrZh, axis=1)>0)/len(partelse_matchrZh), "From matched to ISR, FSR")
        #print(fillne(gentt[g+'pdgId'])[ttpart_matchrZh])
        #######
        self.events  = {'genWeight':self.events['genWeight'][:100]} 
        del self.geninfo, self.fatjets
        # now seperate between sdm, and bbvl and bin in pt
        #pt_bins = [200,300,400,500,700,1000]
        pt_bins = [200,300,450,1000]
        msd_cut  = ((matched_fj['FatJet_msoftdrop'] > 50) & (matched_fj['FatJet_msoftdrop'] < 200))
        bbvl_cut = (matched_fj['FatJet_deepTagMD_bbvsLight'] >= 0.8) 
        msd_num, edges = np.histogram(np.clip(matched_fj['FatJet_pt'][msd_cut],                  0, 750), bins=pt_bins)
        bbvl_num, _    = np.histogram(np.clip(matched_fj['FatJet_pt'][bbvl_cut],                 0, 750), bins=pt_bins)
        tot_num, _     = np.histogram(np.clip(matched_fj['FatJet_pt'][((bbvl_cut) & (msd_cut))], 0, 750), bins=pt_bins)
        den, _         = np.histogram(np.clip(matched_fj['FatJet_pt'],                           0, 750), bins=pt_bins)
        if __name__ == '__main__':
            #self.fatjets = matched_fj[(matched_fj['FatJet_pt']>700)]
            self.fatjets = matched_fj
        self.metaData = {'msd_num':msd_num ,'bbvl_num':bbvl_num,'tot_num':tot_num,'tot_den':den,'pt_bins':pt_bins,'tot_events':sum(den)}

    def get_event_selection(self): # after objects have been defined
        return ( (self.fatjets['FatJet_pt'].counts >= 1) &
                 ((self.jets['Jet_pt'][(self.jets['Jet_btagDeepB'] > cfg.ZHbb_btagWP[self.year])].counts > 2) if 'QCD' not in self.sample else True) &
                 (self.jets['Jet_pt'].counts > 4) &
                 self.get_MET_filter() &
                 (self.pass_goldenjson() == 1 if self.isData else True) &
                 (
                     np.where(self.events['run'] >= 319077, self.get_HEM_veto(), True) if 
                     (self.year == '2018' and self.isData) else True
                 )
             )

    def get_skim(self):
        __out_dict = {'events':pd.DataFrame.from_dict(self.events)}
        if not self.isData:
            __out_dict['metaData'] = self.metaData
        # close root file
        #self.f.close()
        #
        return __out_dict

    def is_a_jet(self):
        return  (
            (self.jets['Jet_pt']       > 30) & 
            (abs(self.jets['Jet_eta']) < 2.4) & 
            ((self.jets['Jet_pt'] > 50) | (self.jets['Jet_puId'] >= 4) ) &
            ( self.jets['Jet_jetId'] >= 2) 
        )

    def is_a_fatjet(self):
        return (
            (self.fatjets['FatJet_pt'] >  200) &
            (abs(self.fatjets['FatJet_eta']) < 2.4) &   
            ( self.fatjets['FatJet_jetId'] >= 2) 
        )
                        

    def set_pre_cuts(self): # to speed things up?
        _c_keys = ['nFatJet','FatJet_pt']
        import concurrent.futures
        from modules.AnaDict import AnaDict
        executor = concurrent.futures.ThreadPoolExecutor()
        _c_vars = AnaDict({k: self.tarray(k, executor=executor) for k in _c_keys})
        _precut = (_c_vars['nFatJet'] >= 1)
        #_precut = ((_c_vars['nFatJet'] >= 1) & ((_c_vars['FatJet_pt']>700).counts > 0))
        return _precut    

    
    def plot_kinem_dist(self):
        import matplotlib.pyplot as plt
        #plt.hist(np.clip(self.fatjets['FatJet_msoftdrop'],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        #fix, ax = plt.subplots()
        #ax.hist(np.clip(self.fatjets['FatJet_msoftdrop'][np.nansum(self.ttpart_matchrZh_0p6, axis=1)>0],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        #ax.set_ylabel('Counts')
        ##plt.xlabel('softdrop mass (Ak8 pt > 700, Gen matched)')
        #ax.set_xlabel('softdrop mass (DR(<0.6,ttsys), Gen matched)')
        #ax.legend()
        ##
        #fix, ax = plt.subplots()
        #ax.hist(np.clip(self.fatjets['FatJet_msoftdrop'][(np.nansum(self.ttpart_matchrZh_0p6, axis=1)==0) & (np.nansum(self.ttpart_matchrZh_1p2, axis=1)>0)],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        #ax.set_ylabel('Counts')
        #ax.set_xlabel('softdrop mass (DR(<1.2 & >0.6,ttsys), Gen matched)')
        #ax.legend()
        #
        fix, ax = plt.subplots()
        ax.hist(np.clip(self.fatjets['FatJet_msoftdrop'][np.nansum(self.ttpart_matchrZh_1p2, axis=1)==0],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        ax.set_ylabel('Counts')
        ax.set_xlabel('softdrop mass (DR(>1.2,ttsys), Gen matched)')
        ax.legend()
        #
        pt_cut = 450
        #fix, ax = plt.subplots()
        #ax.hist(np.clip(self.fatjets['FatJet_msoftdrop'][(self.fatjets['FatJet_pt']>pt_cut) & (np.nansum(self.ttpart_matchrZh_0p6, axis=1)>0)],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        #ax.set_ylabel('Counts')
        ##plt.xlabel('softdrop mass (Ak8 pt > 700, Gen matched)')
        #ax.set_xlabel(f'softdrop mass (DR(<0.6,ttsys), pt>{pt_cut}, Gen matched)')
        #ax.legend()
        ##
        #fix, ax = plt.subplots()
        #ax.hist(np.clip(self.fatjets['FatJet_msoftdrop'][(self.fatjets['FatJet_pt']>pt_cut) & (np.nansum(self.ttpart_matchrZh_0p6, axis=1)==0) & (np.nansum(self.ttpart_matchrZh_1p2, axis=1)>0)],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        #ax.set_ylabel('Counts')
        #ax.set_xlabel(f'softdrop mass (DR(<1.2 & >0.6,ttsys), pt>{pt_cut}, Gen matched)')
        #ax.legend()
        #
        fix, ax = plt.subplots()
        print(self.fatjets['FatJet_pt'][(self.fatjets['FatJet_pt']>pt_cut) & (np.nansum(self.ttpart_matchrZh_1p2, axis=1)==0)])
        ax.hist(np.clip(self.fatjets['FatJet_msoftdrop'][(self.fatjets['FatJet_pt']>pt_cut) & (np.nansum(self.ttpart_matchrZh_1p2, axis=1)==0)],0,249), range=(0,250), histtype='step', label=self.sample+'_'+self.year)
        ax.set_ylabel('Counts')
        ax.set_xlabel(f'softdrop mass (DR(>1.2,ttsys), pt>{pt_cut}, Gen matched)')
        ax.legend()



if __name__ == '__main__':
    year = '2016'
    golden_json=json.load(open(cfg.goodLumis_file[year]))
    #test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed//2017/QCD_HT_2000toInf_2017/302603A9-C719-8247-8F99-4084F124BB01_Skim_5.root'
    #_ = bbvLEffSkim(test_file, 'QCD_HT_2000toInf', year, isData=False, is4eff=True, jec_sys=None, golden_json=golden_json)
    test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed//2016/TTZToBB_2016/F92669C5-93D9-0C48-A241-1DC0F7602E3E_Skim_1.root'
    _ = bbvLEffSkim(test_file, 'TTZToBB', year, isData=False, is4eff=True, jec_sys=None, golden_json=golden_json)
    #test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed//2016/ttHTobb_2016/F74205CA-7EDA-7845-B4DB-470E3D9A93CC_Skim_3.root'
    #_ = bbvLEffSkim(test_file, 'ttHTobb', year, isData=False, is4eff=True, jec_sys=None, golden_json=golden_json)

    print(_.get_skim())
    save_pdf("bbvl_sdm_genmatch_testing.pdf")(_.plot_kinem_dist)()
