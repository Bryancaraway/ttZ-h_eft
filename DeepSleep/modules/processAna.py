########################
### process code     ###
### for analysis     ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#

import time
from numba import njit, prange
from uproot_methods import TLorentzVectorArray
# Custom cfg, lib, modules
from cfg import deepsleepcfg as cfg
import lib.fun_library as lib
from lib.fun_library import fill1e, fillne, deltaR, deltaPhi, invM, calc_mtb
from modules.AnaDict import AnaDict
#
import numpy as np
np.random.seed(0)
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
    b_wp   = cfg.ZHbb_btagWP
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
        self.process_data()

    def process_data(self):
        # to contain lepCleaned_v2, match_gen, and ZHbbAna
        start = time.perf_counter()
        self.lepCleaned_v2()
        # might be better to reconstruct zh here, code kinda messy :(
        if self.isSignal or self.isttbar:
            self.match_gen_lep() 
            if   self.isttbar:
                self.match_gen_tt()  
            elif self.isSignal:
                # Note: can be improved if ran after ZH is reconstructed and used to gen-reco match 
                self.match_gen_sig() 
        #
        self.recoZh()
        #
        self.ak4_df.to_pickle(f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_ak4.pkl")
        self.ak8_df.to_pickle(f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_ak8.pkl")
        self.val_df.to_pickle(f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_val.pkl")
        self.gen_df.to_pickle(f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_gen.pkl")
        self.rtc_df.to_pickle(f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_rtc.pkl")
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
        gen_ids = self.gen_df['GenPart_pdgId']
        gen_mom = self.gen_df['GenPart_genPartIdxMother']
        is_tt_bb = (((abs(gen_ids) == 5) & (abs(gen_ids[gen_mom]) > 6)).sum() >= 2)
        self.val_df['tt_bb'] = is_tt_bb
        
    def match_gen_sig(self):
        # match tt or ttZ/h gen particles to recontructed objects
        g = 'GenPart_'
        gen_ids, gen_mom, gen_pt, gen_eta, gen_phi = self.gen_df.loc(
            [g+'pdgId',g+'genPartIdxMother',g+'pt',g+'eta',g+'phi']
        ).values()
        #
        fj= 'FatJet_'
        ak8_pt, ak8_eta, ak8_phi, ak8_Zhbbtag = self.ak8_df.loc(
            [fj+'pt'+self.lc,fj+'eta'+self.lc,fj+'phi'+self.lc, fj+'btagHbb'+self.lc]
        )[self.ak8_df[fj+'lep_mask']].values()
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
        zh_match_dR = deltaR(zh_eta,zh_phi,ak8_eta,ak8_phi)
        zh_match = ((ak8_pt >= self.pt_cut ) & (ak8_Zhbbtag >= 0.0) &  
                    (zh_match_dR <= 0.8) & (zh_pt >= ( self.pt_cut-100.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
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
        self.val_df['matchedGenLep'] = ((lep_match_dr < .1).sum() > 0)

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
        ak4_pt, ak4_eta, ak4_phi, ak4_E, ak4_btag = self.ak4_df.loc(
            [j+str_+self.lc for str_ in ['pt','eta','phi','E','btagDeepB']]
        )[self.ak4_df[j+'lep_mask']].values()
        b_pt, b_eta, b_phi, b_E, b_btag = [ak4_k[ak4_btag >= self.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_E,ak4_btag]]
        q_pt, q_eta, q_phi, q_E, q_btag = [ak4_k[ak4_btag <  self.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_E,ak4_btag]]
        ak8_pt, ak8_eta, ak8_phi, ak8_E, sd_M, ak8_bbtag, ak8_Zhbbtag, w_tag, t_tag, subj1, subj2 = self.ak8_df.loc(
            [fj+str_+self.lc for str_ in ['pt','eta','phi','E','msoftdrop','btagDeepB','btagHbb','deepTag_WvsQCD','deepTag_TvsQCD','subJetIdx1','subJetIdx2']]
        )[self.ak8_df[fj+'lep_mask']].values()
        subj_pt, subj_btag = self.ak8_df['SubJet_pt'], self.ak8_df['SubJet_btagDeepB']
        #
        lep_pt, lep_eta, lep_phi, lep_E = [self.val_df[key] for key in [l+'pt',l+'eta',l+'phi',l+'E']]
        met_pt, met_phi = self.val_df['MET_pt'], self.val_df['MET_phi']
        # reconstruct z, h -> bb candidate
        Zh_reco_cut = ((ak8_pt >= self.pt_cut) & (sd_M >= 50) & (sd_M <= 200) & (ak8_Zhbbtag >= 0.0)) # in future, put kinem cut values in cfg file
        Zh_Zhbbtag,Zh_pt,Zh_eta,Zh_phi,Zh_E,Zh_M,Zh_wtag,Zh_ttag,Zh_bbtag=lib.sortbyscore(
            [ak8_Zhbbtag,ak8_pt,ak8_eta,ak8_phi,ak8_E,sd_M,w_tag,t_tag,ak8_bbtag],ak8_Zhbbtag,Zh_reco_cut)
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
        Zh_ak4_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],ak4_eta,ak4_phi)
        Zh_b_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],b_eta,b_phi)
        Zh_q_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],q_eta,q_phi)
        Zh_l_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],lep_eta,lep_phi)
        Zh_l_invM_sd = lib.invM_sdM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],lep_pt,lep_eta,lep_phi,lep_E)
        l_b_dr = deltaR(lep_eta.values,lep_phi.values,*map(fillne,[b_eta,b_phi]))
        l_b_invM_E = lib.invM_E(lep_pt.values,lep_eta.values,lep_phi.values,lep_E.values,*map(fillne,[b_pt,b_eta,b_phi,b_E]))
        getTLV  = TLorentzVectorArray.from_ptetaphi
        b_tlv   = getTLV( *map(lambda b : fillne(b).T, [b_pt,b_eta,b_phi,b_E]))
        lep_tlv = getTLV( lep_pt,lep_eta,lep_phi,lep_E)
        bl_tlv = lep_tlv + b_tlv
        lb_mtb = lib.calc_mtb(bl_tlv.pt.T,bl_tlv.phi.T,met_pt.values,met_phi.values)            
        #
        ind_lb = np.argsort(l_b_dr,axis=1) 
        l_b_invM_E_dRsort, l_b_mtb_dRsort, l_b_dr_dRsort = [np.take_along_axis(lb_comb,ind_lb,axis=1) for lb_comb in [l_b_invM_E,lb_mtb,l_b_dr]]
        max_l_b_dr,  min_l_b_dr  = np.nanmax(l_b_dr_dRsort, axis=1),     np.nanmin(l_b_dr_dRsort, axis=1)
        max_lb_invM, min_lb_invM = np.nanmax(l_b_invM_E_dRsort, axis=1), np.nanmin(l_b_invM_E_dRsort, axis=1)
        #
        n_b_outZhbb = (Zh_b_dr > .8).counts
        n_q_outZhbb = (Zh_q_dr > .8).counts
        n_b_inZhbb  = (Zh_b_dr < .8).counts
        n_q_inZhbb  = (Zh_q_dr < .8).counts
        Zh_ak4_dr, Zh_b_dr, Zh_q_dr = fillne(Zh_ak4_dr), fillne(Zh_b_dr), fillne(Zh_q_dr)
        ind_Zh_b = np.argsort(np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis=1)
        b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_E_dRsort, b_btag_dRsort  = [np.take_along_axis(fillne(b_k),ind_Zh_b,axis=1) for b_k in [b_pt,b_eta,b_phi,b_E, b_btag]]
        b1_tlv_dRsort = getTLV(b_pt_dRsort[:,0],b_eta_dRsort[:,0],b_phi_dRsort[:,0],b_E_dRsort[:,0])
        b2_tlv_dRsort = getTLV(b_pt_dRsort[:,1],b_eta_dRsort[:,1],b_phi_dRsort[:,1],b_E_dRsort[:,1])
        b12_pt_dRsort =  (b1_tlv_dRsort+b2_tlv_dRsort).pt
        Zh_b_invM_sd = lib.invM_sdM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0],Zh_M[:,0],b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_E_dRsort)
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
        self.val_df['Zh_E']      = Zh_E[:,0]
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
        self.val_df['l_b1_invM'] = l_b_invM_E_dRsort[:,0]
        self.val_df['l_b1_dr']   = l_b_dr_dRsort[:,0]
        self.val_df['l_b2_mtb']  = l_b_mtb_dRsort[:,1]
        self.val_df['l_b2_invM'] = l_b_invM_E_dRsort[:,1]
        self.val_df['l_b2_dr']   = l_b_dr_dRsort[:,1]
        #####################################################

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

if __name__ == '__main__':
    #
    # will need to open pkl files for testing
    sample = 'DY'
    print('Reading Files...')
    dir_ = 'test_ja_files/'
    ak4_df = AnaDict.read_pickle(dir_+f'{sample}_2017_ak4.pkl')
    ak8_df = AnaDict.read_pickle(dir_+f'{sample}_2017_ak8.pkl')
    val_df = pd     .read_pickle(dir_+f'{sample}_2017_val.pkl')
    gen_df = AnaDict.read_pickle(dir_+f'{sample}_2017_gen.pkl')
    rtc_df = AnaDict.read_pickle(dir_+f'{sample}_2017_rtc.pkl')
    print('Processing data...')
    process_ana_dict = {'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df, 'sample':sample}
    processAna(process_ana_dict)
