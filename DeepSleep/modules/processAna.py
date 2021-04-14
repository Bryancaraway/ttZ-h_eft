########################
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
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import time
import re
import json
from numba import njit, prange
from uproot_methods import TLorentzVectorArray
# Custom cfg, lib, modules
import config.ana_cff as cfg
from config.sample_cff import sample_cfg
import lib.fun_library as lib
from lib.fun_library import fill1e, fillne, deltaR, deltaPhi, invM, calc_mtb, t2Run
from modules.AnaDict import AnaDict
#
import numpy as np
np.random.seed(0)
np.seterr(invalid='ignore')
import pandas as pd
##

eps = 0.0001

class processAna :
    '''
    callable object designed to process data that has been extracted from
    root files
    
    the output are .pkl files ready for histogram analysis and Combine Tool
    '''
    # Allowed class variables
    outDir   = cfg.master_file_path
    dataDir  = cfg.dataDir
    outTag   = ''
    year     = '2017'
    files    = None
    sample   = None
    isData   = False
    isSignal = False
    isttbar  = False
    #
    condor = False
    keep_all = False
    #
    #lc     = '_drLeptonCleaned'
    lc = ''
    pt_cut = cfg.ZHptcut
    b_wp   = cfg.ZHbb_btagWP[year]
    #
    ak4_df  = None
    ak8_df  = None
    gen_df  = None
    meta_df = None
    val_df  = None
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
        #self.lepCleaned_v2()
        #
        #self.recoZh()
        from modules.zh_helper import reco_zh_helper_andrew as recoZH_andrew
        recoZH_andrew(self)
        # dum fix for files w/o this variable
        self.val_df['topptWeight']      = 1.
        self.val_df['topptWeight_Up']   = 1.
        self.val_df['topptWeight_Down'] = 1.
        #
        if self.isSignal or self.isttbar:
            self.match_gen_lep() 
            if   self.isttbar:
                self.match_gen_tt()  
            elif self.isSignal:
                # Note: can be improved if ran after ZH is reconstructed and used to gen-reco match 
                self.match_gen_sig() 
        #
        self.passHLT_by_year()
        #self.applyDNN(model_file='withdak8md_model.h5')
        #self.applyDNN(model_file='noak8md_model.h5')
        self.applyDNN(model_file='withbbvl_model.h5')
        self.applyDNN(model_file='withbbvlnewgenm_model.h5')
        self.val_df['NN'] = self.val_df['withbbvlnewgenm_NN'] # hard-coded to not break anythin

        # add hlt variables into MC and set to True for convience later
        if not self.isData:
            self.finalize_btag_w()
            self.getbbvlsf()
            self.getLepEffSF()
            self.getLepSF()
            #self.getBCbtagSF()
        self.categorize_process()
        #
        if self.outTag != '': outTag = f'{self.outTag}_'
        else: outTag = '' 
        #
        if self.condor: out_path = f'{self.sample}_{self.year}_{outTag}'
        else:
            out_path=f"{self.outDir}{self.year}/{'mc_files' if not self.isData else 'data_files'}/{self.sample}_{outTag}"
        self.val_df.to_pickle(out_path+"val.pkl")
        if self.keep_all:
            self.ak4_df.to_pickle(out_path+"ak4.pkl")
            self.ak8_df.to_pickle(out_path+"ak8.pkl")
            if not self.isData:
                self.gen_df.to_pickle(out_path+"gen.pkl")
            #self.rtc_df.to_pickle(out_path+"rtc.pkl")
        finish = time.perf_counter()
        print(f'\nTime to finish {type(self).__name__} for {self.sample}: {finish-start:.1f}\n')

    def finalize_btag_w(self):
        r_ratio_json = json.load(open(cfg.dataDir+f'/btagw_r_ratio/btagw_r_ratio_{self.year}.json', 'r'))
        if 'sys' in sample_cfg[self.sample]['out_name']:
            # ttbar/ttbb_sys to use nominal ttbar and ttbb
            r_ratios = r_ratio_json[sample_cfg['_'.join(self.sample.split('_')[:-1])]['out_name']] 
        elif 'EFT' in self.sample: # dont really care about b-tag weight here, just take some dummy wieght
            if 'TTZ' in self.sample:
                r_ratios = r_ratio_json['ttZ'] # signal eft
            if 'TTH' in self.sample:
                r_ratios = r_ratio_json['ttH'] # signal eft
            else: 
                r_ratios = r_ratio_json['ttbb'] # ttbb (perhaps ttbar) eft
        else:
            r_ratios = r_ratio_json[sample_cfg[self.sample]['out_name']]
        nj  = self.val_df['n_ak4jets'].clip(0,9)
        for r_key  in r_ratios:
            btw_name = r_key.replace('r_ratio','BTagWeight')
            if btw_name in self.val_df:
                r_ = nj.apply((lambda i : r_ratios[r_key][int(i)]))
                self.val_df[btw_name+'_nocorr'] = self.val_df[btw_name]
                self.val_df.loc[:,btw_name] = r_*self.val_df[btw_name]

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
        # first get ttbar decay type : Had, Semi, Di
        if   '2L2Nu' in self.sample:
            self.val_df['tt_type'] = 'Di'
        elif 'Semi' in self.sample:
            self.val_df['tt_type'] = 'Semi'
        elif 'Had' in self.sample:
            self.val_df['tt_type'] = 'Had'
        # match tt or ttZ/h gen particles to recontructed objects
        g = 'GenPart_' 
        gen_ids, gen_mom, gen_st, gen_pt, gen_eta, gen_phi, gen_mass, gentt_bb = self.gen_df.loc(
            [g+'pdgId',g+'genPartIdxMother',g+'status',g+'pt',g+'eta',g+'phi', g+'mass', 'genTtbarId']
        ).values()
        rZh_eta = self.val_df['Zh_eta'].values
        rZh_phi = self.val_df['Zh_phi'].values
        #
        gentt_bb = gentt_bb % 100
        is_tt_C = ( (gentt_bb>=41) & (gentt_bb<50) )
        is_tt_B = gentt_bb>=51
        print(gentt_bb[is_tt_B])

        is_tt_b  = gentt_bb == 51
        is_tt_2b = gentt_bb == 52
        is_tt_bb = gentt_bb >= 53
        # 
        self.val_df['tt_C'] = is_tt_C
        self.val_df['tt_B'] = is_tt_B
        self.val_df['tt_b'] = is_tt_b
        self.val_df['tt_2b'] = is_tt_2b
        self.val_df['tt_bb'] = is_tt_bb
        print('Total',len(is_tt_B))
        print('tt+B', sum(is_tt_B))
        print('tt+b', sum(is_tt_b))
        print('tt+2b', sum(is_tt_2b))
        print('tt+bb', sum(is_tt_bb))
        # calc invM of extra bb
        #if 'bb' in self.sample.lower():
        ext_bb = (lambda c: c[(gen_mom <= 0) & (abs(gen_ids) == 5)])#[((gen_mom <= 0) & (abs(gen_ids) == 5)).sum() == 2])
        getTLVm = TLorentzVectorArray.from_ptetaphim
        b1 = getTLVm(*map((lambda b: b.pad(2)[:,0]), list(map(ext_bb,[gen_pt,gen_eta,gen_phi, gen_mass]))))
        b2 = getTLVm(*map((lambda b: b.pad(2)[:,1]), list(map(ext_bb,[gen_pt,gen_eta,gen_phi, gen_mass]))))
        self.val_df['ttbb_genbb_invm'] = (b1+b2).mass
        self.val_df['ttbb_genbb_pt'] = (b1+b2).pt
        #
        # calculate toppt weight for powheg only
        tt_pt = gen_pt[(abs(gen_ids) == 6)]
        # Using the newer theo (NNLO QCD + NLO EW) corrections which is better for BSM analysis aspects
        sf = (lambda x: 0.103*np.exp(-0.0118*np.clip(x,0,np.inf)) - 0.000134*np.clip(x,0,np.inf) + 0.973) #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Case_3_3_The_Effective_Field_The

        #toppt_sys = AnaDict.read_pickle(f'{self.dataDir}toppt_sys_files/toppt_sys.pkl')
        #https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf'
        # the theo toppt event re-weighting unc. is based on [1, w**2] where w is the event reweighting 
        toppt_rwgt = np.sqrt(sf(tt_pt[:,0]) * sf(tt_pt[:,1])) 
        toppt_rwgt_up = np.where(toppt_rwgt > 1.0, toppt_rwgt**2,  1.0)
        toppt_rwgt_dn = np.where(toppt_rwgt < 1.0, toppt_rwgt**2,  1.0)
        self.val_df['topptWeight']      = toppt_rwgt
        self.val_df['topptWeight_Up']   = toppt_rwgt_up
        self.val_df['topptWeight_Down'] = toppt_rwgt_dn
            
                    
    @staticmethod
    def getToppt_sys(dist, sys):
        bins=[float(bin_str.split(',')[0]) for bin_str in sys] + [float(list(sys.keys())[-1].split(',')[1])]
        digi = pd.cut(np.clip(dist,bins[0]+eps,bins[-1]-eps), bins=bins, labels= sorted(sys.keys(), key= lambda z : float(z.split(',')[0])))
        top_sys = digi.map(lambda k: sys[k])
        return top_sys.to_numpy(dtype=float)

            
        
    def match_gen_sig(self):
        # match tt or ttZ/h gen particles to recontructed objects
        g = 'GenPart_' 
        gen_ids, gen_mom, gen_st, gen_pt, gen_eta, gen_phi, gen_mass = self.gen_df.loc(
            [g+'pdgId',g+'genPartIdxMother',g+'status',g+'pt',g+'eta',g+'phi',g+'mass']
        ).values()
        #
        #fj= 'FatJet_'
        #ak8_pt, ak8_eta, ak8_phi, ak8_Zhbbtag = self.ak8_df.loc(
        #    [fj+'pt'+self.lc,fj+'eta'+self.lc,fj+'phi'+self.lc, fj+'btagHbb'+self.lc]
        #)[self.ak8_df[fj+'lep_mask']].values()
        rZh_eta = self.val_df['Zh_eta'].values
        rZh_phi = self.val_df['Zh_phi'].values
        #
        isbb_fromZ     = ((abs(gen_ids) == 5) & (gen_ids[gen_mom] == 23))
        isqq_fromZ     = ((abs(gen_ids) <  5) & (gen_ids[gen_mom] == 23))
        isllnunu_fromZ = ((abs(gen_ids) >=  11) & (abs(gen_ids) <= 16) & (gen_ids[gen_mom] == 23))
        isbb_fromH    = ((abs(gen_ids) == 5) & (gen_ids[gen_mom] == 25))
        isnonbb_fromH = ((abs(gen_ids) != 5) & (gen_ids[gen_mom] == 25))
        isHbb     = ((gen_ids == 25) & (isbb_fromH.sum() == 2))
        isHnonbb  = ((gen_ids == 25) & (isbb_fromH.sum() == 0))
        isZbb  = ((gen_ids == 23) & (isbb_fromZ.sum() == 2) & ((isHbb.sum() == 0) & (isHnonbb.sum() == 0)))
        isZqq  = ((gen_ids == 23) & (isqq_fromZ.sum() == 2) & ((isHbb.sum() == 0) & (isHnonbb.sum() == 0)))
        isZllnunu = ((gen_ids == 23) & (isllnunu_fromZ.sum() == 2) & ((isHbb.sum() == 0) & (isHnonbb.sum() == 0)))
        isZH = ((isHbb) | (isZbb) | (isZqq) | (isZllnunu) | (isHnonbb))
        print("isHbb", sum(isHbb.sum()))
        print("isHnonbb", sum(isHnonbb.sum()))
        print("isZbb", sum(isZbb.sum()))
        print("isZqq", sum(isZqq.sum()))
        print("isZllnunu", sum(isZllnunu.sum()))
        print("recoZH",len(rZh_eta),len(rZh_phi))
        #
        zh_pt  = fill1e(gen_pt [isZH]).flatten()
        zh_eta = fill1e(gen_eta[isZH]).flatten()
        zh_phi = fill1e(gen_phi[isZH]).flatten()
        zh_mass = fill1e(gen_mass[isZH]).flatten()
        #
        zh_st  = fill1e(gen_st [isZH]).flatten()
        #print(np.unique(zh_st,return_counts=True))
        #
        zh_match_dR = deltaR(zh_eta,zh_phi,rZh_eta, rZh_phi)
        rzh_matchb_dR = deltaR(rZh_eta,rZh_phi,gen_eta[(isbb_fromZ) | (isbb_fromH)], gen_phi[(isbb_fromZ) | (isbb_fromH)])
        zh_matchbb    = ((rzh_matchb_dR <= 0.6).sum() == 2)
        zh_matchb    = ((rzh_matchb_dR <= 0.6).sum() == 1)
        zh_nomatchb    = ((rzh_matchb_dR <= 0.6).sum() == 0)
        #zh_match_dR = deltaR(zh_eta,zh_phi,ak8_eta,ak8_phi)
        #zh_match = ((ak8_pt >= self.pt_cut ) & (ak8_Zhbbtag >= 0.0) &  
        #            (zh_match_dR <= 0.8) & (zh_pt >= ( self.pt_cut-100.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
        zh_match = ((zh_match_dR <= 0.8) & (zh_pt >= ( self.pt_cut-100.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
        #
        self.val_df['Zbb']= (isZbb.sum() > 0)
        self.val_df['Hbb']= (isHbb.sum() > 0)
        self.val_df['Hnonbb']= (isHnonbb.sum() > 0)
        self.val_df['Zqq']= (isZqq.sum() > 0)
        self.val_df['Zllnunu']= (isZllnunu.sum() > 0)
        self.val_df['genZHpt']  = zh_pt
        self.val_df['genZHeta'] = zh_eta
        self.val_df['genZHphi'] = zh_phi
        self.val_df['genZHmass'] = zh_mass
        getTLVm = TLorentzVectorArray.from_ptetaphim
        zh_tlv   = getTLVm(zh_pt,zh_eta,zh_phi,zh_mass)
        self.val_df['genZHrap'] = zh_tlv.rapidity
        self.val_df['genZHstxs'] = abs(zh_tlv.rapidity) <= 2.5
        #
        self.val_df['matchedGenZH']    = (zh_match).sum() > 0 
        self.val_df['matchedGen_Zbb']  = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isZbb.sum() >  0))
        self.val_df['matchedGen_Hbb']  = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isHbb.sum() >  0))
        self.val_df['matchedGen_ZHbb'] = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & ((isZbb.sum() + isHbb.sum()) > 0))
        self.val_df['matchedGen_Zqq']  = (((zh_match).sum() > 0) & (self.val_df['matchedGenLep']) & (isZqq.sum() >  0))
        #
        self.val_df['matchedGen_ZHbb_bb']  = ((self.val_df['matchedGen_ZHbb'] == 1) & (zh_matchbb  == 1))
        self.val_df['matchedGen_ZHbb_b']   = ((self.val_df['matchedGen_ZHbb'] == 1) & (zh_matchb   == 1))
        self.val_df['matchedGen_ZHbb_nob'] = ((self.val_df['matchedGen_ZHbb'] == 1) & (zh_nomatchb == 1))
        #
        self.val_df['matchedGen_ZHbb_nn']  = ((self.val_df['matchedGen_ZHbb'] == 1) & ((zh_matchbb  == 1)|(zh_matchb   == 1)))
        #
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
        #print('Number of leptons in sample (from W and W from Top) that pass our SingleLepton Req.')
        #print(np.unique(gen_ids[islep].counts[self.val_df['Hbb']], return_counts=True))
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
        #ak4_pt, ak4_eta, ak4_phi, ak4_mass, ak4_btag = self.ak4_df.loc(
        #    [j+str_+self.lc for str_ in ['pt','eta','phi','mass','btagDeepB']]
        #)[self.ak4_df[j+'lep_mask']].values()
        ak4_pt, ak4_eta, ak4_phi, ak4_mass, ak4_btag = self.ak4_df.loc(
            [j+str_+self.lc for str_ in ['pt','eta','phi','mass','btagDeepB']]
        ).values()
        ht= ak4_pt.sum()
        #print(ht)
        self.val_df['HT'] = ht
        b_pt, b_eta, b_phi, b_mass, b_btag = [ak4_k[ak4_btag >= self.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_mass,ak4_btag]]
        q_pt, q_eta, q_phi, q_mass, q_btag = [ak4_k[ak4_btag <  self.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_mass,ak4_btag]]
        #ak8_pt, ak8_eta, ak8_phi, sd_M, ak8_bbtag, ak8_Zhbbtag, w_tag, t_tag, subj1, subj2 = self.ak8_df.loc(
        #    [fj+str_+self.lc for str_ in ['pt','eta','phi','msoftdrop','btagDeepB','deepTagMD_bbvsLight','deepTagMD_WvsQCD','deepTagMD_TvsQCD','subJetIdx1','subJetIdx2']]
        #)[self.ak8_df[fj+'lep_mask']].values()
        ak8_pt, ak8_eta, ak8_phi, sd_M, ak8_bbtag, ak8_Zhbbtag, w_tag, t_tag, subj1, subj2 = self.ak8_df.loc(
            [fj+str_+self.lc for str_ in ['pt','eta','phi','msoftdrop','btagDeepB','deepTagMD_bbvsLight','deepTagMD_WvsQCD','deepTagMD_TvsQCD','subJetIdx1','subJetIdx2']]
        ).values()
        subj_pt, subj_btag = self.ak8_df['SubJet_pt'], self.ak8_df['SubJet_btagDeepB']
        #
        lep_pt, lep_eta, lep_phi, lep_M = [self.val_df[key] for key in [l+'pt',l+'eta',l+'phi',l+'mass']]
        met_pt, met_phi = self.val_df['MET_pt'], self.val_df['MET_phi']
        # reconstruct z, h -> bb candidate
        Zh_reco_cut = ((ak8_pt >= self.pt_cut) & (sd_M >= 50) & (sd_M <= 200) & (ak8_Zhbbtag >= 0.0)) # in future, put kinem cut values in cfg file
        #print(Zh_reco_cut)
        Zh_Zhbbtag,Zh_pt,Zh_eta,Zh_phi,Zh_M,Zh_wtag,Zh_ttag,Zh_bbtag=lib.sortbyscore(
            [ak8_Zhbbtag,ak8_pt,ak8_eta,ak8_phi,sd_M,w_tag,t_tag,ak8_bbtag],ak8_Zhbbtag,Zh_reco_cut)
        # =============================== # 
        # compute subject info related to Zh candidate : have to sort out candidate with no subjet
        ak8_sj1_pt, ak8_sj2_pt, ak8_sj1_btag, ak8_sj2_btag = [subj_k[id_[id_ != -1]] for subj_k in [subj_pt, subj_btag] for id_ in [subj1, subj2]] # sj kinems per ak8jet
        ak8_sj1_sdM, ak8_sj2_sdM, ak8_sj1_Zhbb, ak8_sj2_Zhbb = [ak8_k[id_ != -1] for ak8_k in [sd_M, ak8_Zhbbtag] for id_ in [subj1, subj2]] # ak8 kinems 
        #Zh_kinem_sj1_cut = ((ak8_sj1_pt>=self.pt_cut) & (ak8_sj1_sdM >= 50) & (ak8_sj1_sdM <= 200) & (ak8_sj1_Zhbb >= 0.0))
        Zh_kinem_sj1_cut = (ak8_sj1_Zhbb >= 0.0)
        Zh_sj1_pt, Zh_sj1_btag, Zh_sj1_Zhbb = lib.sortbyscore([
            ak8_sj1_pt,ak8_sj1_btag,ak8_sj1_Zhbb],ak8_sj1_Zhbb,Zh_kinem_sj1_cut)
        #Zh_kinem_sj2_cut = ((ak8_sj2_pt>=self.pt_cut) & (ak8_sj2_sdM >= 50) & (ak8_sj2_sdM <= 200) & (ak8_sj2_Zhbb >= 0.0))
        Zh_kinem_sj2_cut = ((ak8_sj2_Zhbb >= 0.0))
        Zh_sj2_pt, Zh_sj2_btag, Zh_sj2_Zhbb = lib.sortbyscore([
            ak8_sj2_pt,ak8_sj2_btag,ak8_sj2_Zhbb],ak8_sj2_Zhbb,Zh_kinem_sj2_cut)
        #

        Zh_sj_b12  = np.column_stack([Zh_sj1_btag[:,0], Zh_sj2_btag[:,0]])
        Zh_sj_pt12 = np.nan_to_num(np.column_stack([Zh_sj1_pt[:,0], Zh_sj2_pt[:,0]]) )

        Zh_sjpt12_over_fjpt = (Zh_sj_pt12[:,0] +  Zh_sj_pt12[:,1])/Zh_pt[:,0]
        Zh_sjpt1_over_fjpt, Zh_sjpt2_over_fjpt  = [(Zh_sj_pt12[:,i])/Zh_pt[:,0] for i in range(2)]
        # =============================== #
        # Calculate sphericity and aplanarity of the event
        spher, aplan = lib.calc_SandA(
            np.append(fillne(ak4_pt),  lep_pt .to_numpy()[:,np.newaxis], axis=1),
            np.append(fillne(ak4_eta), lep_eta.to_numpy()[:,np.newaxis], axis=1),
            np.append(fillne(ak4_phi), lep_phi.to_numpy()[:,np.newaxis], axis=1)
        )
        # =============================== #
        # Caclulate combinatorix between Zh and ak4, b, q, l || l and b
        Zh_ak4_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(ak4_eta),fillne(ak4_phi))
        Zh_b_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(b_eta),fillne(b_phi))
        ind_Zh_b = np.argsort(np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis=1)
        #ind_Zh_b = np.argsort(-1*np.where(Zh_b_dr > 0.8, fillne(b_btag), np.nan),axis=1)
        b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort, b_btag_dRsort  = [
            np.take_along_axis(np.where(Zh_b_dr > 0.8, fillne(b_k),np.nan),ind_Zh_b,axis=1) for b_k in [b_pt,b_eta,b_phi,b_mass, b_btag]]
        #
        Zh_q_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(q_eta),fillne(q_phi))
        #ind_Zh_q = np.argsort(np.where(Zh_q_dr > 0.8, Zh_q_dr, np.nan),axis=1)
        ind_Zh_q = np.argsort(-1*np.where(Zh_q_dr > 0.8, fillne(q_pt), np.nan),axis=1)
        q_pt_dRsort, q_eta_dRsort, q_phi_dRsort, q_mass_dRsort, q_btag_dRsort  = [
            np.take_along_axis(np.where(Zh_q_dr > 0.8, fillne(q_k),np.nan),ind_Zh_q,axis=1) for q_k in [q_pt,q_eta,q_phi,q_mass, q_btag]]

        Zh_l_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],lep_eta,lep_phi)
        #
        Zh_l_invM_sd = lib.invM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],lep_pt,lep_eta,lep_phi,lep_M)
        #l_b_dr = deltaR(lep_eta.values,lep_phi.values,*map(fillne,[b_eta,b_phi]))
        b_b_dr = deltaR(b_eta_dRsort[:,0],b_phi_dRsort[:,0], b_eta_dRsort[:,1],b_phi_dRsort[:,1])
        q_q_dr = deltaR(q_eta_dRsort[:,0],q_phi_dRsort[:,0], q_eta_dRsort[:,1],q_phi_dRsort[:,1])
        l_b_dr = deltaR(lep_eta.values,lep_phi.values, b_eta_dRsort,b_phi_dRsort)
        #l_b_invM = lib.invM(lep_pt.values,lep_eta.values,lep_phi.values,lep_M.values,*map(fillne,[b_pt,b_eta,b_phi,b_mass]))
        l_b_invM = lib.invM(lep_pt.values,lep_eta.values,lep_phi.values,lep_M.values, b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_mass_dRsort)
        #
        getTLV  = TLorentzVectorArray.from_ptetaphi
        getTLVm = TLorentzVectorArray.from_ptetaphim
        #
        b_tlv   = getTLVm(*map(lambda p: p.T,[b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_mass_dRsort]) )
        q_tlv   = getTLVm(*map(lambda p: p.T,[q_pt_dRsort,q_eta_dRsort,q_phi_dRsort,q_mass_dRsort]) )
        lep_tlv = getTLV( lep_pt.values,lep_eta.values,lep_phi.values,lep_M.values)
        bl_tlv = lep_tlv + b_tlv
        self.val_df['max_lb_invM_alt'] = np.nanmax(bl_tlv.mass.T, axis=1)
        lb_mtb = lib.calc_mtb(bl_tlv.pt.T,bl_tlv.phi.T,met_pt.values,met_phi.values)            
        #
        ind_lb = np.argsort(l_b_dr,axis=1) 
        l_b_invM_dRsort, l_b_mtb_dRsort, l_b_dr_dRsort = [np.take_along_axis(lb_comb,ind_lb,axis=1) for lb_comb in [l_b_invM,lb_mtb,l_b_dr]]
        nearl_b_pt_dRsort, nearl_b_eta_dRsort, nearl_b_phi_dRsort, nearl_b_mass_dRsort = [np.take_along_axis(lb_k,ind_lb,axis=1)[:,0] for lb_k in [b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_mass_dRsort]]
        max_l_b_dr,  min_l_b_dr  = np.nanmax(l_b_dr_dRsort, axis=1),     np.nanmin(l_b_dr_dRsort, axis=1)
        max_lb_invM, min_lb_invM = np.nanmax(l_b_invM_dRsort, axis=1), np.nanmin(l_b_invM_dRsort, axis=1)
        # farthest b to l
        ind_lb_far = np.argsort(-l_b_dr,axis=1) 
        farl_b_pt_dRsort, farl_b_eta_dRsort, farl_b_phi_dRsort, farl_b_mass_dRsort = [np.take_along_axis(lb_k,ind_lb_far,axis=1)[:,0] for lb_k in [b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_mass_dRsort]]
        #print(l_b_dr)
        #print(ind_lb_far)
        far_l_b_q_dr = deltaR(farl_b_eta_dRsort,farl_b_phi_dRsort, q_eta_dRsort,q_phi_dRsort) # get 1st and 2nd closest quarks
        max_far_l_b_q_dr = np.nanmax(far_l_b_q_dr, axis=1) # new
        min_far_l_b_q_dr = np.nanmin(far_l_b_q_dr, axis=1) # new
        ind_farl_bq_dr = np.argsort(far_l_b_q_dr,axis=1)
        far_l_b_q_pt_dRsort, far_l_b_q_eta_dRsort, far_l_b_q_phi_dRsort, far_l_b_q_mass_dRsort = [np.take_along_axis(q_,ind_farl_bq_dr,axis=1) for q_ in [q_pt_dRsort,q_eta_dRsort,q_phi_dRsort,q_mass_dRsort]]
        far_l_b_tlv = getTLVm(farl_b_pt_dRsort,farl_b_eta_dRsort,farl_b_phi_dRsort,farl_b_mass_dRsort)
        near_l_b_tlv = getTLVm(nearl_b_pt_dRsort,nearl_b_eta_dRsort,nearl_b_phi_dRsort,nearl_b_mass_dRsort)
        q1_tlv = getTLVm(far_l_b_q_pt_dRsort[:,0],far_l_b_q_eta_dRsort[:,0],far_l_b_q_phi_dRsort[:,0],far_l_b_q_mass_dRsort[:,0])
        q2_tlv = getTLVm(far_l_b_q_pt_dRsort[:,1],far_l_b_q_eta_dRsort[:,1],far_l_b_q_phi_dRsort[:,1],far_l_b_q_mass_dRsort[:,1])
        bqq_tlv = far_l_b_tlv + q1_tlv + q2_tlv
        bqq_mass = bqq_tlv.mass # new
        Zh_bqq_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],bqq_tlv.eta, bqq_tlv.phi) # new
        lbbqq_tlv = bqq_tlv + lep_tlv + near_l_b_tlv
        Zh_lbbqq_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],lbbqq_tlv.eta, lbbqq_tlv.phi) # new
        #
        n_b_outZhbb = np.nansum(Zh_b_dr > .8, axis=1)
        n_q_outZhbb = np.nansum(Zh_q_dr > .8 , axis=1)
        n_b_inZhbb  = np.nansum(Zh_b_dr <= .8, axis=1)
        n_q_inZhbb  = np.nansum(Zh_q_dr <= .8, axis=1)

        #Zh_ak4_dr, Zh_b_dr, Zh_q_dr = fillne(Zh_ak4_dr), fillne(Zh_b_dr), fillne(Zh_q_dr)

        b1_tlv_dRsort = getTLVm(b_pt_dRsort[:,0],b_eta_dRsort[:,0],b_phi_dRsort[:,0],b_mass_dRsort[:,0])
        b2_tlv_dRsort = getTLVm(b_pt_dRsort[:,1],b_eta_dRsort[:,1],b_phi_dRsort[:,1],b_mass_dRsort[:,1])
        b12_pt_dRsort =  (b1_tlv_dRsort+b2_tlv_dRsort).pt
        b12_m_dRsort  =  (b1_tlv_dRsort+b2_tlv_dRsort).mass # new
        b12_dr_dRsort =  b1_tlv_dRsort.delta_r(b2_tlv_dRsort) # new 
        h_t_b         = np.nansum(b_pt_dRsort, axis=1) # new
        #
        sc_pt_outZh   = h_t_b + np.nansum(q_pt_dRsort, axis=1) + lep_pt # new
        #
        Zh_b_invM_sd = lib.invM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0],Zh_M[:,0],b_pt_dRsort,b_eta_dRsort,b_phi_dRsort,b_mass_dRsort) 
        mtb1 = calc_mtb(b_pt_dRsort[:,0],b_phi_dRsort[:,0],met_pt,met_phi)
        mtb2 = calc_mtb(b_pt_dRsort[:,1],b_phi_dRsort[:,1],met_pt,met_phi)
        best_Wb_invM_sd = np.where(((mtb2 > mtb1) & (mtb2 != np.nan)), Zh_b_invM_sd[:,1], Zh_b_invM_sd[:,0]) 
        # find best resolved top candidate from ak4 jets outside of Zh candidate
        ak4_outZh= np.where(Zh_ak4_dr>=.8,Zh_ak4_dr,np.nan)
        #
        self.val_df['nBottoms']  = b_pt.counts
        self.val_df['n_ak4jets'] = ak4_pt.counts
        self.val_df['n_ak8jets'] = ak8_pt.counts
        #
        self.val_df['jetpt_1']  = ak4_pt.pad(2)[:,0]
        self.val_df['jetpt_2']  = ak4_pt.pad(2)[:,1]
        self.val_df['bjetpt_1']  = b_pt.pad(2)[:,0]
        self.val_df['bjetpt_2']  = b_pt.pad(2)[:,1]
        self.val_df['jeteta_1'] = ak4_eta.pad(2)[:,0]
        self.val_df['jeteta_2'] = ak4_eta.pad(2)[:,1]
        self.val_df['bjeteta_1'] = b_eta.pad(2)[:,0]
        self.val_df['bjeteta_2'] = b_eta.pad(2)[:,1]
        self.val_df['jetbtag_1'] = ak4_btag.pad(2)[:,0]
        self.val_df['jetbtag_2'] = ak4_btag.pad(2)[:,1]
        self.val_df['bjetbtag_1'] = b_btag.pad(2)[:,0]
        self.val_df['bjetbtag_2'] = b_btag.pad(2)[:,1]
        #
        self.val_df['fjetpt_1']  = ak8_pt.pad(1)[:,0]
        self.val_df['fjeteta_1'] = ak8_eta.pad(1)[:,0]
        self.val_df['fjetsdm_1'] = sd_M.pad(1)[:,0]
        self.val_df['fjetbbvl_1']= ak8_Zhbbtag.pad(1)[:,0]
        self.val_df['fjetwscore_1']= w_tag.pad(1)[:,0]
        self.val_df['fjettscore_1']= t_tag.pad(1)[:,0]
        #
        self.val_df['n_ak8_Zhbb'] = ak8_pt[((ak8_pt >= self.pt_cut) & (sd_M >= 50) & (sd_M <= 200) & (ak8_Zhbbtag >= 0.8))].counts
        self.val_df['Zh_bbvLscore']  = Zh_Zhbbtag[:,0]
        self.val_df['Zh_pt']     = Zh_pt[:,0]
        self.val_df['Zh_eta']    = Zh_eta[:,0]
        self.val_df['Zh_phi']    = Zh_phi[:,0]
        self.val_df['Zh_M']      = Zh_M[:,0]
        self.val_df['Zh_Wscore'] = Zh_wtag[:,0]
        self.val_df['Zh_Tscore'] = Zh_ttag[:,0]
        self.val_df['outZh_max_Wscore'] = np.max(np.nan_to_num(Zh_wtag[:,1:]), axis=1) # new
        self.val_df['outZh_max_Tscore'] = np.max(np.nan_to_num(Zh_ttag[:,1:]), axis=1) # new
        self.val_df['outZh_max_bbvLscore']  = np.max(np.nan_to_num(Zh_Zhbbtag[:,1:]), axis=1) # new
        self.val_df['Zh_deepB']  = Zh_bbtag[:,0]
        self.val_df['n_Zh_btag_sj'] = np.nansum(Zh_sj_b12 >= self.b_wp, axis=1)
        self.val_df['n_Zh_sj']       = np.nansum(Zh_sj_b12 >= 0, axis=1)
        self.val_df['Zh_bestb_sj']   = np.nan_to_num(np.nanmax(Zh_sj_b12, axis=1))
        self.val_df['Zh_worstb_sj']  = np.nan_to_num(np.nanmin(Zh_sj_b12, axis=1))
        self.val_df['Zh_bbscore_sj'] = np.nansum(np.nan_to_num(Zh_sj_b12), axis=1)
        self.val_df['sjpt12_over_Zhpt'] = Zh_sjpt12_over_fjpt
        self.val_df['sjpt1_over_Zhpt'] = Zh_sjpt1_over_fjpt
        self.val_df['sjpt2_over_Zhpt'] = Zh_sjpt2_over_fjpt

        self.val_df['spher'] = spher
        self.val_df['aplan'] = aplan
        
        self.val_df['outZh_bb_dr'] = b_b_dr #new
        self.val_df['outZh_qq_dr'] = np.nan_to_num(q_q_dr) #new     

        self.val_df['max_lb_dr'] = max_l_b_dr
        self.val_df['min_lb_dr'] = min_l_b_dr
        self.val_df['max_lb_invM'] = max_lb_invM
        self.val_df['min_lb_invM'] = min_lb_invM

        self.val_df['b1_outZh_score'] = b_btag_dRsort[:,0]
        self.val_df['b2_outZh_score'] = b_btag_dRsort[:,1]
        self.val_df['b1_over_Zhpt']       = b_pt_dRsort[:,0]/Zh_pt[:,0]
        self.val_df['b2_over_Zhpt']       = b_pt_dRsort[:,1]/Zh_pt[:,0]
        self.val_df['bb_over_Zhpt']       = b12_pt_dRsort/Zh_pt[:,0]
        self.val_df['b12_outZh_m'] = b12_m_dRsort
        self.val_df['b12_outZh_dr'] = b12_dr_dRsort
        self.val_df['ht_b'] = h_t_b # new
        self.val_df['ht_outZh'] = sc_pt_outZh # new
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
        self.val_df['l_b1_mtb']  = l_b_mtb_dRsort[:,0]
        self.val_df['l_b1_invM'] = l_b_invM_dRsort[:,0]
        self.val_df['l_b1_dr']   = l_b_dr_dRsort[:,0]
        self.val_df['l_b2_mtb']  = l_b_mtb_dRsort[:,1]
        self.val_df['l_b2_invM'] = l_b_invM_dRsort[:,1]
        self.val_df['l_b2_dr']   = l_b_dr_dRsort[:,1]
        #
        self.val_df['max_farl_b_q_dr'] = np.nan_to_num(max_far_l_b_q_dr) # new
        self.val_df['min_farl_b_q_dr'] = np.nan_to_num(min_far_l_b_q_dr) # new
        self.val_df['outZh_bqq_mass'] = np.nan_to_num(bqq_mass) # new
        self.val_df['Zh_bqq_dr'] = np.nan_to_num(Zh_bqq_dr) # new
        self.val_df['Zh_lbbqq_dr'] = np.nan_to_num(Zh_lbbqq_dr) # new
        #
    
    def applyDNN(self, model_file='withdak8md_model.h5'):
        from modules.dnn_model import DNN_model

        nn_dir    = cfg.dnn_ZH_dir 
        dnn_vars = {
            #'withdak8md_model.h5' :cfg.allvars_dnn_ZH_vars,
            #'noak8md_model.h5'     :cfg.nodak8md_dnn_ZH_vars,
            #'selvars_ttzh_model.h5':cfg.selvars_dnn_ZH_vars,
            'withbbvl_model.h5':    cfg.withbbvl_dnn_ZH_vars,
            #'ttzh_newgenm.h5':      cfg.withbbvl_dnn_ZHgenm_vars,
            'withbbvlnewgenm_model.h5':      cfg.withbbvl_dnn_ZHgenm_vars,
        }
        resetIndex = (lambda df: df.reset_index(drop=True).copy())
        #open json nn settings
        # --
        #m_info =  {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 'other_settings': {'fl_a': [0.75, 1, 1.25], 'fl_g': 0.25, 'lr_alpha': 0.0003}, 'n_epochs': 80, 'batch_size': 10256}
        #m_info = {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 'other_settings': {'fl_a': [1, 0.75, 0.75], 'fl_g': 0.25, 'lr_alpha': 0.0003}, 'n_epochs': 100, 'batch_size': 10256}
        m_info =  {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 'other_settings': {'fl_a': [0.75, 1, 0.25], 'fl_g': 0.25, 'lr_alpha': 0.0003}, 'n_epochs': 100, 'batch_size': 10256}
        #
        dnn = DNN_model(m_info['sequence'],m_info['other_settings']) 
        #nn_model = dnn.Build_Model(len(cfg.dnn_ZH_vars), load_weights='ttzh_model.h5')#'nn_ttzh_model.h5')
        nn_model = dnn.Build_Model(len(dnn_vars[model_file]), load_weights=model_file)#'nn_ttzh_model.h5')
        #nn_model = dnn.Build_Model(len(cfg.dnn_ZH_vars), load_weights='nn_ttzh_model.h5')
        #
        base_cuts = lib.getZhbbBaseCuts(self.val_df)
        #
        pred_df = resetIndex(self.val_df[dnn_vars[model_file]][base_cuts])

        pred = nn_model.predict(pred_df)[:,2]
        #print(pred)
        df_nn_name = model_file.replace("model.h5",'')+'NN'
        self.val_df.loc[:,df_nn_name] = -1.
        self.val_df.loc[base_cuts,df_nn_name] = pred

        #
    def addHLT_to_MC(self):
        # add hlt variables into MC and set them to 1
        for var in cfg.ana_vars['dataHLT_all']+cfg.ana_vars[f'dataHLT_{self.year}']:
            self.val_df[var] = True

    def getbbvlsf(self):
        sf_file_dict = json.load(open(cfg.dataDir+f'/deepak8sf/deepak8_bbvL_sf_{self.year}.json','r'))
        pt_bins = sf_file_dict['pt_bins']
        score_bins = sf_file_dict['score_bins']
        # apply by file 
        def get_right_sf(_df,r_or_f='fake'):
            sf = sf_file_dict[r_or_f]['score_pt_sf']
            sf_err = sf_file_dict[r_or_f]['score_pt_sf_err']
            pt_digi = pd.cut(_df['Zh_pt'].clip(pt_bins[0]+0.001,pt_bins[-1]-0.001), bins = pt_bins, 
                             right=True,  include_lowest=True, labels=range(len(pt_bins)-1))
            score_digi = pd.cut(_df['Zh_bbvLscore'], bins = score_bins, 
                                right=True,  include_lowest=True, labels=range(len(score_bins)-1))
            pt_digi_ignore_nan    = np.nan_to_num(pt_digi.values).astype(int)
            score_digi_ignore_nan = np.nan_to_num(score_digi.values).astype(int)
            
            get_sf = (lambda _sf : np.array([ _sf[x][y] for x,y in zip(score_digi_ignore_nan,pt_digi_ignore_nan)]))
            _df.loc[:,'dak8md_bbvl_sf'] = np.where(score_digi.isna().values, 1, get_sf(sf))     # for bbvl = -10, we drop these events later
            list_sf_err                 = np.where(score_digi.isna().values, 0, get_sf(sf_err)) # for bbvl = -10
            _df.loc[:,'dak8md_bbvl_sf_Up']   = _df['dak8md_bbvl_sf'] + list_sf_err
            _df.loc[:,'dak8md_bbvl_sf_Down'] = _df['dak8md_bbvl_sf'] - list_sf_err
            return _df
            
        if self.isSignal:
            self.val_df.loc[:,'dak8md_bbvl_sf']      = 1
            self.val_df.loc[:,'dak8md_bbvl_sf_Up']   = 1
            self.val_df.loc[:,'dak8md_bbvl_sf_Down'] = 1
            self.val_df = pd.concat([get_right_sf(self.val_df.loc[((self.val_df['Hbb']==True)|(self.val_df['Zbb']==True))], 'real'),
                                     get_right_sf(self.val_df.loc[((self.val_df['Hbb']==False)&(self.val_df['Zbb']==False))],'fake')
                                 ],axis='rows')
        elif 'TTbb' in self.sample:
            self.val_df = get_right_sf(self.val_df,'real')
        else:
            self.val_df = get_right_sf(self.val_df,'fake')
            

    @t2Run
    def getLepEffSF(self):
        self.val_df['lep_trigeffsf']      = 1
        self.val_df['lep_trigeffsf_Up']   = 1
        self.val_df['lep_trigeffsf_Down'] = 1
        for lep, singlelep in zip(['Electron','Muon'],['passSingleLepElec','passSingleLepMu']):
            lep_pt, lep_eta  = self.val_df['Lep_pt'], self.val_df['Lep_eta']
            if lep == 'Muon': lep_eta = abs(lep_eta)
            sf_dict = json.load(open(cfg.dataDir+'/lep_effsf_files/'+f"trigeffSF_{lep}_{self.year}.json",'r'))[lep]
            sf, sf_up, sf_down, pt_bins, eta_bins = sf_dict['pt_eta_sf'], sf_dict['pt_eta_sf_Up'], sf_dict['pt_eta_sf_Down'], sf_dict['pt_bins'], sf_dict['eta_bins']
            pt_digi  = pd.cut(lep_pt.clip(pt_bins[0]+eps,pt_bins[-1]-eps), bins=pt_bins, right=True, include_lowest=True, labels=range(len(pt_bins)-1))
            eta_digi = pd.cut(np.clip(lep_eta,eta_bins[0]+eps,eta_bins[-1]-eps), bins=eta_bins, right=True, include_lowest=True, labels=range(len(eta_bins)-1))
            get_sf = (lambda _sf : np.array([ _sf[x][y] for x,y in zip(pt_digi.values,eta_digi.values)]))
            weight, weight_up, weight_down = get_sf(sf), get_sf(sf_up), get_sf(sf_down)
            self.val_df.loc[self.val_df[singlelep]==True, 'lep_trigeffsf']      = (weight)            [self.val_df[singlelep] == True]
            self.val_df.loc[self.val_df[singlelep]==True, 'lep_trigeffsf_Up']   = (weight+weight_up)  [self.val_df[singlelep] == True]
            self.val_df.loc[self.val_df[singlelep]==True, 'lep_trigeffsf_Down'] = (weight-weight_down)[self.val_df[singlelep] == True]
            #print(self.val_df['lep_trigeffsf'])
            #print(self.val_df['lep_trigeffsf_Up'])
            #print(self.val_df['lep_trigeffsf_Down'])

    @t2Run
    def getLepSF(self):
        lep_pt  = self.val_df['Lep_pt']
        lep_eta = self.val_df['Lep_eta'] 
        for lep_type in ('muon','electron'):
            #sf = AnaDict.read_pickle(f'{self.outDir}{self.year}/lepton_sf_files/{lep_type}_sf_{self.year}.pkl')
            sf = AnaDict.read_pickle(f'{self.dataDir}lep_sf_files/{lep_type}_sf_{self.year}.pkl')
            total_lep_sf = [np.ones(len(lep_pt)),np.zeros(len(lep_pt))]
            #total_lep_sf = [np.ones(len(lep_pt)),np.ones(len(lep_pt))]
            for k in sf:
                # dict structure is sf type : eta_bins: pt_bins: values, up, down
                if ((lep_type == 'muon' and self.year != '2016') or k == 'SF'): lep_eta = abs(lep_eta)
                eta_keys = list(sf[k].keys())
                pt_keys  = list(sf[k][eta_keys[-1]].keys())
                eta_bins=[float(bin_str.split(',')[0]) for bin_str in eta_keys] + [float(eta_keys[-1].split(',')[1])]
                pt_bins =[float(bin_str.split(',')[0]) for bin_str in pt_keys] +  [float( pt_keys[-1].split(',')[1])]
                eta_digi = pd.cut(np.clip(lep_eta,eta_bins[0]+eps,eta_bins[-1]-eps), bins=eta_bins, right=True, include_lowest=True,
                                  labels=sorted(eta_keys, key= lambda z : float(z.split(',')[0])))
                pt_digi  = pd.cut(lep_pt.clip(pt_bins[0]+eps,pt_bins[-1]-eps), bins=pt_bins, right=True, include_lowest=True, 
                                  labels=sorted(pt_keys, key= lambda z : float(z.split(',')[0])))
                #lep_sf = pd.concat([eta_digi,pt_digi], axis='columns').apply(lambda x: sf[k][x[0]][x[1]]['values'], axis='columns')
                lep_sf = np.array([(lambda x,y: np.array([sf[k][x][y]['values'], sf[k][x][y]['up']]))(x,y) for x,y in zip(eta_digi.values,pt_digi.values)]) # this should be propagating the error properly
                #total_lep_sf *= lep_sf.to_numpy(dtype=float)
                total_lep_sf[0] *= lep_sf[:,0]
                lep_sf_err = abs(lep_sf[:,1] - lep_sf[:,0]) * lep_sf[:,0]
                total_lep_sf[1] = total_lep_sf[0]*np.sqrt(np.power(lep_sf_err/lep_sf[:,0],2)+ np.power(total_lep_sf[1]/total_lep_sf[0],2)) 
            self.val_df[f'{lep_type}_sf']      = total_lep_sf[0]
            if lep_type == 'muon': 
                total_lep_sf[1] = np.sqrt(total_lep_sf[1]**2 + (.03*total_lep_sf[0])**2) # 3% additional error added to suslep mu
            self.val_df[f'{lep_type}_sf_up']   = total_lep_sf[1]+total_lep_sf[0]
            self.val_df[f'{lep_type}_sf_down'] = total_lep_sf[0]-total_lep_sf[1]
        #
        for var in ['','_up','_down']:
            self.val_df[f'lep_sf{var}'] = pd.concat(
                [self.val_df[f'electron_sf{var}'][self.val_df['passSingleLepElec']==1],
                 self.val_df[f'muon_sf{var}'][self.val_df['passSingleLepMu']==1]]).sort_index()

    def addLepEff(self, dist, lep_type, eff):
        bins=[float(bin_str.split(',')[0]) for bin_str in eff] + [float(list(eff.keys())[-1].split(',')[1])]
        digi = pd.cut(dist.clip(bins[0]+eps,bins[-1]-eps), bins=bins, labels= sorted(eff.keys(), key= lambda z : float(z.split(',')[0])))
        lep_eff = digi.map(lambda k: eff[k]['values'])
        lep_eff_up = digi.map(lambda k: eff[k]['up'])
        lep_eff_down = digi.map(lambda k: eff[k]['down'])
        return (lep_eff.to_numpy(dtype=float),lep_eff_up.to_numpy(dtype=float),lep_eff_down.to_numpy(dtype=float))
     
    #@t2Run
    #def getBCbtagSF(self):
    #    sf = cfg.BC_btag_sf[self.year]['values']
    #    n_b = self.val_df['nBottoms_drLeptonCleaned'].clip(0,6)
    #    BCbtagSF = n_b.map(lambda k: sf[int(k)])
    #    self.val_df['BC_btagSF'] = BCbtagSF
        
    def passHLT_by_year(self):
        self.val_df['pbt_elec'] = cfg.hlt_path['electron'][self.year](self.val_df)
        self.val_df['pbt_muon'] = cfg.hlt_path['muon'][self.year](self.val_df)
        self.val_df['isEleE']  = (self.val_df['pbt_elec'] & self.val_df['passSingleLepElec'])
        self.val_df['isMuonE'] = (self.val_df['pbt_muon'] & self.val_df['passSingleLepMu'])
    
    def categorize_process(self):
        # categorize_processes in sample
        # for datacard making script
        #
        # special rules for handling TTZH, TTZ_bb
        # TTBar, tt_bb, and Data
        self.val_df['process'] = ''
        
        def handleTTZ():
            self.val_df['process'] = 'ttZ'
            if self.sample == 'TTZToQQ':
                self.val_df.loc[(self.val_df['Zbb']== True),'process'] = 'old_ttZbb'
        def handleTTH():
            self.val_df['process'] = 'ttH'
        def handleTTBar():
            self.val_df.loc[(self.val_df['tt_B'] != True),'process'] = 'TTBar'
            self.val_df.loc[(self.val_df['tt_B'] == True),'process'] = 'old_tt_B'
            self.add_tt_C_rate_unc()
        def handleTT_bb():
            self.val_df.loc[(self.val_df['tt_B'] != True),'process'] = 'non_tt_B'
            self.val_df.loc[(self.val_df['tt_B'] == True),'process'] = 'tt_B' 
            #self.val_df.loc[(self.val_df['tt_2b']== True),'process'] = 'tt_2b' # this is a subset of tt_B 
            self.add_weights_to_ttbb()
            self.add_tt_2b_rate_unc()
        def handleST():
            self.val_df['process'] = 'single_t'
        def handleTTX():
            self.val_df['process'] = 'ttX'
        def handleVjets():
            self.val_df['process'] = 'VJets'
        def handleOther():
            self.val_df['process'] = 'other'
        def handleEleData():
            self.val_df.loc[((self.val_df['pbt_elec'] == True) & 
                             (self.val_df['passSingleLepElec'] == True)),'process'] = 'Data'
            self.val_df.loc[(self.val_df['process'] != 'Data'), 'process'] = 'non_Data'
        def handleMuData():
            self.val_df.loc[((self.val_df['pbt_muon'] == True) & 
                             (self.val_df['passSingleLepMu'] == True)),'process'] =   'Data'
            self.val_df.loc[(self.val_df['process'] != 'Data'), 'process'] = 'non_Data'

        sample_to_process = {'ttZ'                                  : handleTTZ,
                             'ttH'                                  : handleTTH,
                             'TTZ_EFT'                              : handleTTZ,
                             'TTH_EFT'                              : handleTTH,
                             'TTBar'                                : handleTTBar,
                             'TTJets_EFT'                           : handleTTBar,
                             'ttbb'                                 : handleTT_bb,
                             'TTbb_EFT'                             : handleTT_bb,
                             'single_t'                             : handleST,
                             'ttX'                                  : handleTTX,
                             'VJets'                                : handleVjets,
                             'VV'                                   : handleOther,
                             'VVV'                                  : handleOther,
                             'Data_SingleElectron'                  : handleEleData,
                             'Data_SingleMuon'                      : handleMuData
        }
        #sample_to_process.get(self.sample, handleOther)()
        process = sample_cfg[self.sample]['out_name'].replace('_sys','') # replace to handle TTBar_sys, ttbb_sys
        sample_to_process.get(process, handleOther)()
        print(np.unique(self.val_df['process']))

    def add_weights_to_ttbb(self):
        tt_bb_nw_files = f'{cfg.dataDir}/process_norms/process_norms_ttbbw_run2.json'
        tt_bb_norm_weights =  json.load(open(tt_bb_nw_files,'r'))[self.year]
        ttbb_key = (sample_cfg[self.sample]['out_name'], self.sample)
        try:
            self.val_df['weight'] = tt_bb_norm_weights[ttbb_key[0]][ttbb_key[1]]['weight']
        except KeyError: # defualt
            self.val_df['weight'] = tt_bb_norm_weights['ttbb']['TTbb_SemiLeptonic']['weight']
        #print(self.val_df['weight'])
    def add_tt_2b_rate_unc(self):
        self.val_df['tt2bxsecWeight'] = 1
        self.val_df['tt2bxsecWeight_Up']   = np.where((self.val_df['tt_2b']== True), 1.5, 1)
        self.val_df['tt2bxsecWeight_Down'] = np.where((self.val_df['tt_2b']== True), 0.5, 1)
    def add_tt_C_rate_unc(self):
        self.val_df['ttCxsecWeight'] = 1
        self.val_df['ttCxsecWeight_Up']   = np.where((self.val_df['tt_C']== True), 1.5, 1)
        self.val_df['ttCxsecWeight_Down'] = np.where((self.val_df['tt_C']== True), 0.5, 1)


if __name__ == '__main__':

    from modules.AnaDict import AnaDict
    # will need to open pkl files for testing
    sample = 'TTBarLep_pow'
    #sample = 'TT_bb_pow'
    print('Reading Files...')
    dir_ = 'files/2017/mc_files/'
    ak4_df = AnaDict.read_pickle(dir_+f'{sample}_ak4.pkl')
    ak8_df = AnaDict.read_pickle(dir_+f'{sample}_ak8.pkl')
    val_df = pd     .read_pickle(dir_+f'{sample}_val.pkl')
    gen_df = AnaDict.read_pickle(dir_+f'{sample}_gen.pkl')
    rtc_df = AnaDict.read_pickle(dir_+f'{sample}_rtc.pkl')
    print('Processing data...')
    process_ana_dict = {'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df, 'sample':sample, 'year':'2017', 'isData':False, 'isSignal': False, 'isttbar':True}
    processAna(process_ana_dict)

