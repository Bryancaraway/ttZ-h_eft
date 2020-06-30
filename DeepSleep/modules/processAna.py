########################AAA
### process code     ###
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
    sys.path.insert(1,'/home/bcaraway/ttZh_ana/DeepSleep')
import time
import re
from numba import njit, prange
from uproot_methods import TLorentzVectorArray
# Custom cfg, lib, modules
import cfg.deepsleepcfg  as cfg
import lib.fun_library as lib
from lib.fun_library import fill1e, fillne, deltaR, deltaPhi, invM, calc_mtb, t2Run
#
import numpy as np
np.random.seed(0)
np.seterr(invalid='ignore')
import pandas as pd
##

class processAna :
    '''
    callable object designed to process data that has been extracted from
    root files
    
    the output are .pkl files ready for histogram analysis and Combine Tool
    '''
    # Allowed class variables
    outDir   = 'files/'
    year     = '2017'
    files    = None
    sample   = None
    isData   = False
    isSignal = False
    isttbar  = False
    # 
    lc     = '_drLeptonCleaned'
    pt_cut = cfg.ZHptcut
    b_wp   = cfg.ZHbb_btagWP[year]
    #
    ak4_df = None
    ak8_df = None
    gen_df = None
    val_df = None
    rtc_df = None
    #
    ak4_dR = 0.4
    ak8_dR = 0.8

    def __init__(self, kwargs):
        [setattr(self,k,v) for k,v in kwargs.items() if k in processAna.__dict__.keys()]
        self.b_wp = cfg.ZHbb_btagWP[self.year]
        self.process_data()

    def process_data(self):
        # to contain lepCleaned_v2, match_gen, and ZHbbAna
        start = time.perf_counter()
        self.lepCleaned_v2()
        #
        self.recoZh()
        #
        if self.isSignal or self.isttbar:
            self.match_gen_lep() 
            if   self.isttbar:
                self.match_gen_tt()  
            elif self.isSignal:
                # Note: can be improved if ran after ZH is reconstructed and used to gen-reco match 
                self.match_gen_sig() 
        #
        self.applyDNN()
        # add hlt variables into MC and set to True for convience later
        if not self.isData:
            self.addHLT_to_MC()
        self.passHLT_by_year()
        #
        out_path=f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_"
        self.ak4_df.to_pickle(out_path+"ak4.pkl")
        self.ak8_df.to_pickle(out_path+"ak8.pkl")
        self.val_df.to_pickle(out_path+"val.pkl")
        if not self.isData:
            self.gen_df.to_pickle(out_path+"gen.pkl")
        self.rtc_df.to_pickle(out_path+"rtc.pkl")
        finish = time.perf_counter()
        print(f'\nTime to finish {type(self).__name__} for {self.sample}: {finish-start:.1f}\n')

    def lepCleaned_v2(self):
        # Clean ak4 and ak8 jets within dR of 0.4 and 0.8 of Lepton
        lep_eta = self.val_df['Lep_eta'].to_numpy()
        lep_phi = self.val_df['Lep_phi'].to_numpy()
        iter_   = zip([self.ak4_df,self.ak8_df],['Jet','FatJet'], [self.ak4_dR, self.ak8_dR])
        for j_df, j_name, eta_cut in iter_:
            j_eta, j_phi = j_df[j_name+'_eta'+self.lc], j_df[j_name+'_phi'+self.lc]
            lep_j_dR = deltaR(lep_eta,lep_phi, j_eta, j_phi)
            j_df[f'{j_name}_lep_mask'] = (lep_j_dR > eta_cut) # for the future figure out how to apply this and overwrite relevent jet variables
        #
    def match_gen_tt(self):
        # match tt or ttZ/h gen particles to recontructed objects
        #gen_ids = self.gen_df['GenPart_pdgId']
        #gen_mom = self.gen_df['GenPart_genPartIdxMother']
        #gen_st  = self.gen_df['GenPart_status']
        g = 'GenPart_' 
        gen_ids, gen_mom, gen_st, gen_eta, gen_phi = self.gen_df.loc(
            [g+'pdgId',g+'genPartIdxMother',g+'status',g+'eta',g+'phi']
        ).values()
        rZh_eta = self.val_df['Zh_eta'].values
        rZh_phi = self.val_df['Zh_phi'].values
        #is_tt_bb = (((abs(gen_ids) == 5) & (abs(gen_ids[gen_mom]) > 6)).sum() >= 2)
        extra_bb  = ((abs(gen_ids) == 5) & (abs(gen_ids[gen_mom]) != 6) & ((gen_st == 23) | (gen_st == 21)))
        is_tt_bb  = ((extra_bb).sum() >= 2)
        extra_bb_dr = deltaR(rZh_eta,rZh_phi,gen_eta,gen_phi)
        self.val_df['tt_bb'] = is_tt_bb
        print(np.nansum( self.val_df['tt_bb']))
        self.val_df['genMatched_tt_bb'] = (((extra_bb_dr < 1.2) & (extra_bb)).sum() >= 2)
        print(np.nansum(self.val_df['genMatched_tt_bb']))
        
    def match_gen_sig(self):
        # match tt or ttZ/h gen particles to recontructed objects
        g = 'GenPart_' 
        gen_ids, gen_mom, gen_pt, gen_eta, gen_phi = self.gen_df.loc(
            [g+'pdgId',g+'genPartIdxMother',g+'pt',g+'eta',g+'phi']
        ).values()
        #
        #fj= 'FatJet_'
        #ak8_pt, ak8_eta, ak8_phi, ak8_Zhbbtag = self.ak8_df.loc(
        #    [fj+'pt'+self.lc,fj+'eta'+self.lc,fj+'phi'+self.lc, fj+'btagHbb'+self.lc]
        #)[self.ak8_df[fj+'lep_mask']].values()
        rZh_eta = self.val_df['Zh_eta'].values
        rZh_phi = self.val_df['Zh_phi'].values
        #
        isbb_fromZ = ((abs(gen_ids) == 5) & (gen_ids[gen_mom] == 23))
        isqq_fromZ = ((abs(gen_ids) <  5) & (gen_ids[gen_mom] == 23))
        isbb_fromH = ((abs(gen_ids) == 5) & (gen_ids[gen_mom] == 25))
        isHbb  = ((gen_ids == 25) & (isbb_fromH.sum() == 2))
        isZbb  = ((gen_ids == 23) & (isbb_fromZ.sum() == 2))
        isZqq  = ((gen_ids == 23) & (isqq_fromZ.sum() == 2))
        isZH = ((isHbb) | (isZbb) | (isZqq))
        #
        zh_pt  = fill1e(gen_pt [isZH]).flatten()
        zh_eta = fill1e(gen_eta[isZH]).flatten()
        zh_phi = fill1e(gen_phi[isZH]).flatten()
        #
        zh_match_dR = deltaR(zh_eta,zh_phi,rZh_eta, rZh_phi)
        #zh_match_dR = deltaR(zh_eta,zh_phi,ak8_eta,ak8_phi)
        #zh_match = ((ak8_pt >= self.pt_cut ) & (ak8_Zhbbtag >= 0.0) &  
        #            (zh_match_dR <= 0.8) & (zh_pt >= ( self.pt_cut-100.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
        zh_match = ((zh_match_dR <= 0.8) & (zh_pt >= ( self.pt_cut-100.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
        #
        self.val_df['Zbb']= (isZbb.sum() > 0)
        self.val_df['Hbb']= (isHbb.sum() > 0)
        self.val_df['Zqq']= (isZqq.sum() > 0)
        self.val_df['genZHpt']  = zh_pt
        self.val_df['genZHeta'] = zh_eta
        self.val_df['genZHphi'] = zh_phi
        #
        self.val_df['matchedGenZH']    = (zh_match).sum() > 0 
        self.val_df['matchedGen_Zbb']  = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isZbb.sum() >  0))
        self.val_df['matchedGen_Hbb']  = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isHbb.sum() >  0))
        self.val_df['matchedGen_ZHbb'] = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isZqq.sum() == 0))
        self.val_df['matchedGen_Zqq']  = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isZqq.sum() >  0))

    def match_gen_lep(self):
        lep_eta = self.val_df['Lep_eta'].to_numpy()
        lep_phi = self.val_df['Lep_phi'].to_numpy()
        #
        gen_ids = self.gen_df['GenPart_pdgId']
        gen_mom = self.gen_df['GenPart_genPartIdxMother']
        gen_eta = self.gen_df['GenPart_eta']
        gen_phi = self.gen_df['GenPart_phi']
        #
        islep   = (((abs(gen_ids) == 11) | (abs(gen_ids) == 13)) & ((abs(gen_ids[gen_mom[gen_mom]]) == 6) & (abs(gen_ids[gen_mom]) ==24)))
        lep_match_dr = deltaR(lep_eta,lep_phi,gen_eta[islep],gen_phi[islep])
        self.val_df['matchedGenLep'] = ((lep_match_dr <= .1).sum() > 0)

    def recoZh(self):
        # Note: in order to take advatage of numpy methods with dimensionallity
        # it is necessary to numpy-ize (rectangularize) jagged-arrays using custom (hacky)
        # function 'fillne'
        #
        # create event level variables here
        # get ak4, ak8,val variables needed for analysis #
        l = 'Lep_'
        j = 'Jet_'
        fj= 'FatJet_'
        ak4_pt, ak4_eta, ak4_phi, ak4_mass, ak4_btag = self.ak4_df.loc(
            [j+str_+self.lc for str_ in ['pt','eta','phi','mass','btagDeepB']]
        )[self.ak4_df[j+'lep_mask']].values()
        b_pt, b_eta, b_phi, b_mass, b_btag = [ak4_k[ak4_btag >= self.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_mass,ak4_btag]]
        q_pt, q_eta, q_phi, q_mass, q_btag = [ak4_k[ak4_btag <  self.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_mass,ak4_btag]]
        ak8_pt, ak8_eta, ak8_phi, sd_M, ak8_bbtag, ak8_Zhbbtag, w_tag, t_tag, subj1, subj2 = self.ak8_df.loc(
            [fj+str_+self.lc for str_ in ['pt','eta','phi','msoftdrop','btagDeepB','btagHbb','deepTag_WvsQCD','deepTag_TvsQCD','subJetIdx1','subJetIdx2']]
        )[self.ak8_df[fj+'lep_mask']].values()
        subj_pt, subj_btag = self.ak8_df['SubJet_pt'], self.ak8_df['SubJet_btagDeepB']
        #
        lep_pt, lep_eta, lep_phi, lep_E = [self.val_df[key] for key in [l+'pt',l+'eta',l+'phi',l+'E']]
        met_pt, met_phi = self.val_df['MET_pt'], self.val_df['MET_phi']
        # reconstruct z, h -> bb candidate
        Zh_reco_cut = ((ak8_pt >= self.pt_cut) & (sd_M >= 50) & (sd_M <= 200) & (ak8_Zhbbtag >= 0.0)) # in future, put kinem cut values in cfg file
        Zh_Zhbbtag,Zh_pt,Zh_eta,Zh_phi,Zh_M,Zh_wtag,Zh_ttag,Zh_bbtag=lib.sortbyscore(
            [ak8_Zhbbtag,ak8_pt,ak8_eta,ak8_phi,sd_M,w_tag,t_tag,ak8_bbtag],ak8_Zhbbtag,Zh_reco_cut)
        # compute subject info related to Zh candidate : have to sort out candidate with no subjet
        ak8_sj1_pt, ak8_sj2_pt, ak8_sj1_btag, ak8_sj2_btag = [subj_k[id_[id_ != -1]] for subj_k in [subj_pt, subj_btag] for id_ in [subj1, subj2]]
        ak8_sj1_sdM, ak8_sj2_sdM, ak8_sj1_Zhbb, ak8_sj2_Zhbb = [ak8_k[id_ != -1] for ak8_k in [sd_M, ak8_Zhbbtag] for id_ in [subj1, subj2]]
        Zh_kinem_sj1_cut = ((ak8_sj1_pt>=self.pt_cut) & (ak8_sj1_sdM >= 50) & (ak8_sj1_sdM <= 200) & (ak8_sj1_Zhbb >= 0.0))
        Zh_sj1_pt, Zh_sj1_btag, Zh_sj1_Zhbb = lib.sortbyscore([
            ak8_sj1_pt,ak8_sj1_btag,ak8_sj1_Zhbb],ak8_sj1_Zhbb,Zh_kinem_sj1_cut)
        Zh_kinem_sj2_cut = ((ak8_sj2_pt>=self.pt_cut) & (ak8_sj2_sdM >= 50) & (ak8_sj2_sdM <= 200) & (ak8_sj2_Zhbb >= 0.0))
        Zh_sj2_pt, Zh_sj2_btag, Zh_sj2_Zhbb = lib.sortbyscore([
            ak8_sj2_pt,ak8_sj2_btag,ak8_sj2_Zhbb],ak8_sj2_Zhbb,Zh_kinem_sj2_cut)
        #
        Zh_sj_b12  = np.column_stack([Zh_sj1_btag[:,0], Zh_sj2_btag[:,0]])
        Zh_sj_pt12 = np.nan_to_num(np.column_stack([Zh_sj1_pt[:,0], Zh_sj2_pt[:,0]]) )
        Zh_sjpt12_over_fjpt = (Zh_sj_pt12[:,0] +  Zh_sj_pt12[:,1])/Zh_pt[:,0]
        Zh_sjpt1_over_fjpt, Zh_sjpt2_over_fjpt  = [(Zh_sj_pt12[:,i])/Zh_pt[:,0] for i in range(2)]
        # Calculate sphericity and aplanarity of the event
        spher, aplan = lib.calc_SandA(
            np.append(fillne(ak4_pt),  lep_pt .to_numpy()[:,np.newaxis], axis=1),
            np.append(fillne(ak4_eta), lep_eta.to_numpy()[:,np.newaxis], axis=1),
            np.append(fillne(ak4_phi), lep_phi.to_numpy()[:,np.newaxis], axis=1)
        )
        # Caclulate combinatorix between Zh and ak4, b, q, l || l and b
        Zh_ak4_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(ak4_eta),fillne(ak4_phi))
        Zh_b_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(b_eta),fillne(b_phi))
        Zh_q_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(q_eta),fillne(q_phi))
        Zh_l_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],lep_eta,lep_phi)
        Zh_l_invM_sd = lib.invM_sdM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],lep_pt,lep_eta,lep_phi,lep_E)
        l_b_dr = deltaR(lep_eta.values,lep_phi.values,*map(fillne,[b_eta,b_phi]))
        l_b_invM = lib.invM_Em(lep_pt.values,lep_eta.values,lep_phi.values,lep_E.values,*map(fillne,[b_pt,b_eta,b_phi,b_mass]))
        getTLV  = TLorentzVectorArray.from_ptetaphi
        getTLVm = TLorentzVectorArray.from_ptetaphim
        b_tlv   = getTLVm( *map(lambda b : fillne(b).T, [b_pt,b_eta,b_phi,b_mass])) 
        lep_tlv = getTLV( lep_pt,lep_eta,lep_phi,lep_E)
        bl_tlv = lep_tlv + b_tlv
        lb_mtb = lib.calc_mtb(bl_tlv.pt.T,bl_tlv.phi.T,met_pt.values,met_phi.values)            
        #
        ind_lb = np.argsort(l_b_dr,axis=1) 
        l_b_invM_dRsort, l_b_mtb_dRsort, l_b_dr_dRsort = [np.take_along_axis(lb_comb,ind_lb,axis=1) for lb_comb in [l_b_invM,lb_mtb,l_b_dr]]
        max_l_b_dr,  min_l_b_dr  = np.nanmax(l_b_dr_dRsort, axis=1),     np.nanmin(l_b_dr_dRsort, axis=1)
        max_lb_invM, min_lb_invM = np.nanmax(l_b_invM_dRsort, axis=1), np.nanmin(l_b_invM_dRsort, axis=1)
        #
        n_b_outZhbb = np.nansum(Zh_b_dr > .8, axis=1)
        n_q_outZhbb = np.nansum(Zh_q_dr > .8 , axis=1)
        n_b_inZhbb  = np.nansum(Zh_b_dr <= .8, axis=1)
        n_q_inZhbb  = np.nansum(Zh_q_dr <= .8, axis=1)

        #Zh_ak4_dr, Zh_b_dr, Zh_q_dr = fillne(Zh_ak4_dr), fillne(Zh_b_dr), fillne(Zh_q_dr)

        ind_Zh_b = np.argsort(np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis=1)
        b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort, b_btag_dRsort  = [np.take_along_axis(fillne(b_k),ind_Zh_b,axis=1) for b_k in [b_pt,b_eta,b_phi,b_mass, b_btag]]
        b1_tlv_dRsort = getTLVm(b_pt_dRsort[:,0],b_eta_dRsort[:,0],b_phi_dRsort[:,0],b_mass_dRsort[:,0])
        b2_tlv_dRsort = getTLVm(b_pt_dRsort[:,1],b_eta_dRsort[:,1],b_phi_dRsort[:,1],b_mass_dRsort[:,1])
        b12_pt_dRsort =  (b1_tlv_dRsort+b2_tlv_dRsort).pt
        Zh_b_invM_sd = lib.invM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0],Zh_M[:,0],b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_mass_dRsort) 
        mtb1 = calc_mtb(b_pt_dRsort[:,0],b_phi_dRsort[:,0],met_pt,met_phi)
        mtb2 = calc_mtb(b_pt_dRsort[:,1],b_phi_dRsort[:,1],met_pt,met_phi)
        best_Wb_invM_sd = np.where(((mtb2 > mtb1) & (mtb2 != np.nan)), Zh_b_invM_sd[:,1], Zh_b_invM_sd[:,0]) 
        # find best resolved top candidate from ak4 jets outside of Zh candidate
        ak4_outZh= np.where(Zh_ak4_dr>=.8,Zh_ak4_dr,np.nan)
        rt_disc = fillne(self.rtc_df['ResolvedTopCandidate_discriminator'])
        rt_id1, rt_id2, rt_id3   = self.rtc_df.loc(['ResolvedTopCandidate_j1Idx','ResolvedTopCandidate_j2Idx','ResolvedTopCandidate_j3Idx']).values()
        clean_rtc_fromZh = self.clean_rtc_Zh(rt_disc, rt_id1, rt_id2, rt_id3, ak4_outZh, np.full(rt_disc.shape, 0.0))
        best_rt_score = np.nanmax(clean_rtc_fromZh, axis=1)
        #
        self.val_df['n_ak8_Zhbb'] = ak8_Zhbbtag[Zh_reco_cut].counts
        self.val_df['Zh_score']  = Zh_Zhbbtag[:,0]
        self.val_df['Zh_pt']     = Zh_pt[:,0]
        self.val_df['Zh_eta']    = Zh_eta[:,0]
        self.val_df['Zh_phi']    = Zh_phi[:,0]
        self.val_df['Zh_M']      = Zh_M[:,0]
        self.val_df['Zh_Wscore'] = Zh_wtag[:,0]
        self.val_df['Zh_Tscore'] = Zh_ttag[:,0]
        self.val_df['Zh_deepB']  = Zh_bbtag[:,0]
        self.val_df['n_Zh_btag_sj'] = np.sum(Zh_sj_b12 >= self.b_wp, axis=1)
        self.val_df['n_Zh_sj']       = np.sum(Zh_sj_b12 >= 0, axis=1)
        self.val_df['Zh_bestb_sj']   = np.nan_to_num(np.nanmax(Zh_sj_b12, axis=1))
        self.val_df['Zh_worstb_sj']  = np.nan_to_num(np.min(Zh_sj_b12, axis=1))
        self.val_df['Zh_bbscore_sj'] = np.sum(np.nan_to_num(Zh_sj_b12), axis=1)
        self.val_df['sjpt12_over_Zhpt'] = Zh_sjpt12_over_fjpt
        self.val_df['sjpt1_over_Zhpt'] = Zh_sjpt1_over_fjpt
        self.val_df['sjpt2_over_Zhpt'] = Zh_sjpt2_over_fjpt

        self.val_df['spher'] = spher
        self.val_df['aplan'] = aplan
        
        self.val_df['max_lb_dr'] = max_l_b_dr
        self.val_df['min_lb_dr'] = min_l_b_dr
        self.val_df['max_lb_invM'] = max_lb_invM
        self.val_df['min_lb_invM'] = min_lb_invM

        self.val_df['b1_outZh_score'] = b_btag_dRsort[:,0]
        self.val_df['b2_outZh_score'] = b_btag_dRsort[:,1]
        self.val_df['b1_over_Zhpt']       = b_pt_dRsort[:,0]/Zh_pt[:,0]
        self.val_df['b2_over_Zhpt']       = b_pt_dRsort[:,1]/Zh_pt[:,0]
        self.val_df['bb_over_Zhpt']       = b12_pt_dRsort/Zh_pt[:,0]
        self.val_df['best_Zh_b_invM_sd']  = best_Wb_invM_sd
        self.val_df['Zh_b1_invM_sd'] = Zh_b_invM_sd[:,0]
        self.val_df['Zh_b2_invM_sd'] = Zh_b_invM_sd[:,1]
        self.val_df['nonZhbb_q1_pt'] = -np.sort(-np.where(Zh_q_dr > 0.8, fillne(q_pt), np.nan),axis = 1 )[:,0]
        self.val_df['nonZhbb_q2_pt'] = -np.sort(-np.where(Zh_q_dr > 0.8, fillne(q_pt), np.nan),axis = 1 )[:,1]
        self.val_df['nonZhbb_q1_dr'] =  np.sort( np.where(Zh_q_dr > 0.8, Zh_q_dr, np.nan),axis = 1 )[:,0]
        self.val_df['nonZhbb_q2_dr'] =  np.sort( np.where(Zh_q_dr > 0.8, Zh_q_dr, np.nan),axis = 1 )[:,1]
        self.val_df['nonZhbb_b1_dr'] =  np.sort( np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis = 1 )[:,0]
        self.val_df['nonZhbb_b2_dr'] =  np.sort( np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis = 1 )[:,1]
        self.val_df['n_b_outZh']     = n_b_outZhbb
        self.val_df['n_b_inZh']      = n_b_inZhbb
        self.val_df['n_q_outZh']     = n_q_outZhbb
        self.val_df['n_q_inZh']      = n_q_inZhbb
        self.val_df['Zh_l_dr']       = Zh_l_dr
        self.val_df['Zh_l_invM_sd']  = Zh_l_invM_sd
        self.val_df['best_rt_score'] = best_rt_score
        self.val_df['l_b1_mtb']  = l_b_mtb_dRsort[:,0]
        self.val_df['l_b1_invM'] = l_b_invM_dRsort[:,0]
        self.val_df['l_b1_dr']   = l_b_dr_dRsort[:,0]
        self.val_df['l_b2_mtb']  = l_b_mtb_dRsort[:,1]
        self.val_df['l_b2_invM'] = l_b_invM_dRsort[:,1]
        self.val_df['l_b2_dr']   = l_b_dr_dRsort[:,1]
        #
    
    def applyDNN(self):
        from modules.dnn_model import DNN_model as dnn
        nn_dir   = cfg.dnn_ZH_dir 
        train_df = pd.read_pickle(nn_dir+'train.pkl')
        train_df = pd.read_pickle(nn_dir+'train.pkl')
        trainX   = dnn.resetIndex(train_df.drop(columns=[ 'Signal',*re.findall(r'\w*weight', ' '.join(train_df.keys()))]))
        nn_model = dnn.Build_Model(len(trainX.keys()),1,trainX.mean().values,trainX.std().values)
        #
        base_cuts = lib.getZhbbBaseCuts(self.val_df)
        #
        pred_df = self.val_df[cfg.dnn_ZH_vars][base_cuts]
        pred = nn_model.predict(pred_df.values).flatten()
        
        self.val_df.loc[:,'NN'] = -1.
        self.val_df.loc[base_cuts,'NN'] = pred
        #
    def addHLT_to_MC(self):
        # add hlt variables into MC and set them to 1
        for var in cfg.ana_vars['dataHLT_all']+cfg.ana_vars[f'dataHLT_{self.year}']:
            self.val_df[var] = True

    def passHLT_by_year(self):
        # determine if data/MC passes trigger criteria by year
        elec_dict = {
            '2016': (lambda df: ((df['HLT_Ele27_WPTight_Gsf']==True) | (df['HLT_Photon175']==True) | (df['HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165']==True) | (df['HLT_Ele115_CaloIdVT_GsfTrkIdT']==True))),
            '2017': (lambda df: ((df['HLT_Ele32_WPTight_Gsf_L1DoubleEG']==True) | (df['HLT_Ele35_WPTight_Gsf']==True) | (df['HLT_Photon200']==True))),
            '2018': (lambda df: ((df['HLT_Ele32_WPTight_Gsf']==True) | (df['HLT_Ele115_CaloIdVT_GsfTrkIdT']==True) | (df['HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165']==True) | (df['HLT_Photon200']==True)))
        }
        muon_dict = {
            '2016': (lambda df: ((df['HLT_IsoMu24']==True) | (df['HLT_IsoTkMu24']==True) | (df['HLT_Mu50']==True))),
            '2017': (lambda df: ((df['HLT_IsoMu27']==True) | (df['HLT_Mu50']==True))), #| (df[]==True) | (df[]==True)),
            '2018': (lambda df: ((df['HLT_IsoMu24']==True) | (df['HLT_Mu50']==True) | (df['HLT_OldMu100']==True) | (df['HLT_TkMu100']==True)))
        }
        self.val_df['pbt_elec'] = elec_dict[self.year](self.val_df)
        self.val_df['pbt_muon'] = muon_dict[self.year](self.val_df)

    @staticmethod
    @njit(parallel=True)
    def clean_rtc_Zh(rtd, j1, j2, j3, ak4, out):
        rows, cols = out.shape
        for i in prange(rows):
            for j in prange(cols):
                if np.isnan(j1[i,j]) | np.isnan(j2[i,j]) | np.isnan(j3[i,j]):
                    continue
                if np.isnan(ak4[i,int(j1[i,j])]) | np.isnan(ak4[i,int(j2[i,j])]) | np.isnan(ak4[i,int(j3[i,j])]):
                    continue
                else:
                    out[i,j] = rtd[i,j]
        return out        
        #

if __name__ == '__main__':

    from modules.AnaDict import AnaDict
    # will need to open pkl files for testing
    sample = 'TTBarLep'
    print('Reading Files...')
    dir_ = 'files/2017/mc_files/'
    ak4_df = AnaDict.read_pickle(dir_+f'{sample}_ak4.pkl')
    ak8_df = AnaDict.read_pickle(dir_+f'{sample}_ak8.pkl')
    val_df = pd     .read_pickle(dir_+f'{sample}_val.pkl')
    gen_df = AnaDict.read_pickle(dir_+f'{sample}_gen.pkl')
    rtc_df = AnaDict.read_pickle(dir_+f'{sample}_rtc.pkl')
    print('Processing data...')
    process_ana_dict = {'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df, 'sample':sample, 'year':'2017', 'isData':False, 'isSignal': False, 'isttbar':True, 'outDir': 'files/'}
    processAna(process_ana_dict)
