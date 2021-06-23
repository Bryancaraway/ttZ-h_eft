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
    b_lwp   = cfg.ZHbb_loosebtagWP[year]
    #
    ak4_df    = None
    ak8_df    = None
    gen_df    = None
    meta_df   = None
    val_df    = None
    softe_df  = None
    softmu_df = None
    #
    ak4_dR = 0.4
    ak8_dR = 0.8
    

    def __init__(self, kwargs):
        [setattr(self,k,v) for k,v in kwargs.items() if k in processAna.__dict__.keys()]
        self.b_wp = cfg.ZHbb_btagWP[self.year]
        self.process_data()

    def process_data(self):
        start = time.perf_counter()
        #
        from modules.lep_helper import reco_soft_lep_helper as anti_soft_lep
        anti_soft_lep(self)
        from modules.zh_helper import reco_zh_helper as recoZH
        recoZH(self)
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
        self.applyDNN(model_file='withbbvlnewgenm_model.h5')
        self.applyDNN(model_file='newgenm_model.h5')
        self.applyDNN(model_file='newreduced1p0_model.h5')
        #self.applyDNN(model_file='reduced1p0_model.h5')
        #self.applyDNN(model_file='reduced1p25_model.h5')
        #self.applyDNN(model_file='reduced1p5_model.h5')
        #
        self.val_df['NN'] = self.val_df['withbbvlnewgenm_NN'] # hard-coded to not break anythin
        #self.val_df['NN'] = self.val_df['newgenm_NN'] # hard-coded to not break anythin
        # add hlt variables into MC and set to True for convience later
        if not self.isData:
            self.match_to_rZH()
            self.finalize_btag_w()
            #self.getbbvlsf()
            self.getbbvlsf_v2()
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

    @t2Run
    def match_to_rZH(self):
        g = 'GenPart_'
        gen_ids, gen_mom, gen_eta, gen_phi = self.gen_df.loc(
            [g+'pdgId',g+'genPartIdxMother',g+'eta',g+'phi']
        ).values()
        rZh_eta = self.val_df['Zh_eta'].values
        rZh_phi = self.val_df['Zh_phi'].values
        #
        isb = ( (abs(gen_ids) == 5) & (abs(gen_ids[gen_mom]) != 5) ) 
        b_eta, b_phi =  gen_eta[isb], gen_phi[isb]
        rzh_matchb_dR = deltaR(rZh_eta, rZh_phi, b_eta, b_phi)
        self.val_df['n_genb_matchZH'] = (rzh_matchb_dR <= 0.8).sum()
        self.val_df['bbvl_genmatch'] = ((rzh_matchb_dR <= 0.8).sum() == 2)

    @t2Run
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
    @t2Run
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

            
    @t2Run
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
        istt = (abs(gen_ids) == 6)
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
        rzh_matchtt_dR = deltaR(rZh_eta,rZh_phi,gen_eta[(istt)], gen_phi[(istt)])
        zh_matchbb    = ((rzh_matchb_dR <= 0.8).sum() == 2)
        zh_matchb    = ((rzh_matchb_dR <= 0.8).sum() == 1)
        zh_nomatchb    = ((rzh_matchb_dR <= 0.8).sum() == 0)
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
        self.val_df['matchedGen_ZHbb_tt']  = np.where(self.val_df['matchedGenLep']== True, (rzh_matchtt_dR<=0.8).sum(), -1)
        #
    @t2Run
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
        self.val_df['n_tt_leps'] = (((abs(gen_ids) == 11) | (abs(gen_ids) == 13) | (abs(gen_ids) == 15)) & ((abs(gen_ids[gen_mom[gen_mom]]) == 6) & (abs(gen_ids[gen_mom]) ==24))).sum()
        self.val_df['n_tt_leps_notau'] = (((abs(gen_ids) == 11) | (abs(gen_ids) == 13)) & ((abs(gen_ids[gen_mom[gen_mom]]) == 6) & (abs(gen_ids[gen_mom]) ==24))).sum()
        self.val_df['matchedGenLep'] = ((lep_match_dr <= .1).sum() > 0)
    
    @t2Run
    def applyDNN(self, model_file='withdak8md_model.h5'):
        from modules.dnn_model import DNN_model

        nn_dir    = cfg.dnn_ZH_dir 
        dnn_vars = {
            'newgenm_model.h5' :      cfg.withbbvl_dnn_ZHgenm_vars, # new training
            'newreduced1p0_model.h5' : cfg.reduced1p0genm_vars,
            'withbbvlnewgenm_model.h5':      cfg.withbbvl_dnn_ZHgenm_vars, # 50 var set, old training
            # old
            'reduced1p0_model.h5' : cfg.oldreduced1p0genm_vars,
            'reduced1p25_model.h5' : cfg.oldreduced1p25genm_vars,
            'reduced1p5_model.h5' : cfg.oldreduced1p5genm_vars,
        }
        resetIndex = (lambda df: df.reset_index(drop=True).copy())
        #open json nn settings
        # --
        m_info =  {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 
                   'other_settings': {'fl_a': [0.75, 1, 0.25], 
                                      'fl_g': 0.25, 'lr_alpha': 0.0003}, 
                   'n_epochs': 100, 'batch_size': 10256
        }
        #
        dnn = DNN_model(m_info['sequence'],m_info['other_settings']) 
        nn_model = dnn.Build_Model(len(dnn_vars[model_file]), load_weights=model_file)#'nn_ttzh_model.h5')
        #
        base_cuts = lib.getZhbbBaseCuts(self.val_df)
        #
        pred_df = resetIndex(self.val_df[dnn_vars[model_file]][base_cuts])
        print(self.sample, pred_df)
        if len(pred_df) > 0:
            pred = nn_model.predict(pred_df)[:,2]
        else: 
            pred = 0
        df_nn_name = model_file.replace("model.h5",'')+'NN'
        self.val_df.loc[:,df_nn_name] = -1.
        self.val_df.loc[base_cuts,df_nn_name] = pred
        # only for current main NN
        if model_file == 'withbbvlnewgenm_model.h5': # old, but may keep
            #if model_file == 'newgenm_model.h5': # new
            explain_model = dnn.Build_New_Model(len(dnn_vars[model_file]), nn_model) # output size 64
            if len(pred_df) > 0:
                explain_pred = explain_model.predict(pred_df) # shape (n_events, 64)
            else: 
                explain_pred = -10* np.ones((1,m_info['sequence'][-2][-1]))
            print(m_info['sequence'][-2][-1])
            for i in range(m_info['sequence'][-2][-1]):
                self.val_df.loc[:,f'NN_{i}'] = -10
                self.val_df.loc[base_cuts,f'NN_{i}'] = explain_pred[:,i]
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
    def getbbvlsf_v2(self):
        sf_file_dict = json.load(open(cfg.dataDir+f'/deepak8sf/deepak8_bbvL_sf_{self.year}.json','r'))
        sdm_bins = sf_file_dict['sdm_bins']
        pt_bins = sf_file_dict['pt_bins']
        score_bins = sf_file_dict['score_bins']
        #
        def get_right_sf(_df,r_or_f='fake'):
            out_df = pd.DataFrame()
            sf = sf_file_dict[r_or_f]['score_pt_sdm_sf']
            sf_err = sf_file_dict[r_or_f]['score_pt_sdm_sf_err']
            sdm_digi   = pd.cut(_df['Zh_M'], bins = sdm_bins, 
                                right=True,  include_lowest=True, labels=range(len(sdm_bins)-1))
            pt_digi    = pd.cut(_df['Zh_pt'].clip(pt_bins[0]+0.001,pt_bins[-1]-0.001), bins = pt_bins, 
                                right=True,  include_lowest=True, labels=range(len(pt_bins)-1))
            score_digi = pd.cut(_df['Zh_bbvLscore'], bins = score_bins, 
                                right=True,  include_lowest=True, labels=range(len(score_bins)-1))
            sdm_digi_ignore_nan    = np.nan_to_num(sdm_digi.values).astype(int)
            pt_digi_ignore_nan     = np.nan_to_num(pt_digi.values).astype(int)
            score_digi_ignore_nan  = np.nan_to_num(score_digi.values).astype(int)
            
            get_sf = (lambda _sf : np.array([ _sf[x][y][z] for x,y,z in zip(score_digi_ignore_nan,pt_digi_ignore_nan,sdm_digi_ignore_nan)]))
            out_df['dak8md_bbvl_sf'] = np.where(score_digi.isna().values, 1, get_sf(sf))     # for bbvl = -10, we drop these events later
            list_sf_err                 = np.where(score_digi.isna().values, 0, get_sf(sf_err)) # for bbvl = -10
            out_df['dak8md_bbvl_sf_Up']   = out_df['dak8md_bbvl_sf'] + list_sf_err
            out_df['dak8md_bbvl_sf_Down'] = out_df['dak8md_bbvl_sf'] - list_sf_err
            return out_df

            
        self.val_df.loc[:,['dak8md_bbvl_sf','dak8md_bbvl_sf_Up', 'dak8md_bbvl_sf_Down']] = 1
        self.val_df.loc[self.val_df['bbvl_genmatch'] == True,['dak8md_bbvl_sf','dak8md_bbvl_sf_Up', 'dak8md_bbvl_sf_Down']] = get_right_sf(self.val_df[self.val_df['bbvl_genmatch'] == True], 'real').values
        self.val_df.loc[self.val_df['bbvl_genmatch'] == False,['dak8md_bbvl_sf','dak8md_bbvl_sf_Up', 'dak8md_bbvl_sf_Down']] = get_right_sf(self.val_df[self.val_df['bbvl_genmatch'] == False], 'fake').values


    @t2Run
    def getLepEffSF(self):
        l_dict = {'muon':'Mu', 'electron':'Elec'} 
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
            #
        for l in l_dict:
            for ud in ['','_Up','_Down']:
                self.val_df[f'{l}_trigeffsf{ud}'] = np.where(self.val_df[f'passSingleLep{l_dict[l]}']==1, self.val_df[f'lep_trigeffsf{ud}'], 1)

    @t2Run
    def getLepSF(self):
        l_dict = {'muon':'Mu', 'electron':'Elec'} 
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
            self.val_df[f'{lep_type}_sf']      = np.where(self.val_df[f'passSingleLep{l_dict[lep_type]}']==1, total_lep_sf[0], 1)
            if lep_type == 'muon': 
                total_lep_sf[1] = np.sqrt(total_lep_sf[1]**2 + (.03*total_lep_sf[0])**2) # 3% additional error added to suslep mu
            self.val_df[f'{lep_type}_sf_up']   = np.where(self.val_df[f'passSingleLep{l_dict[lep_type]}']==1, total_lep_sf[1]+total_lep_sf[0], 1)
            self.val_df[f'{lep_type}_sf_down'] = np.where(self.val_df[f'passSingleLep{l_dict[lep_type]}']==1, total_lep_sf[0]-total_lep_sf[1], 1)
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
            #self.val_df['LHE_HTIncoming'] = self.gen_df['LHE_HTIncoming'].flatten()
            #self.val_df['LHE_HT'] = self.gen_df['LHE_HT'].flatten()

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
        def handleRare():
            self.val_df['process'] = 'rare'
        def handleQCD():
            self.val_df['process'] = 'QCD'
            

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
                             'Data_SingleMuon'                      : handleMuData,
                             'rare'                                 : handleRare,
                             'QCD'                                  : handleQCD,
        }
        #sample_to_process.get(self.sample, handleOther)()
        self.val_df["sample"] = self.sample
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

