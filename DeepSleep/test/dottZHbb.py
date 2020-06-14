 #                      #  
##                    ##
########################                               
### Process TTZ/H,   ###
### Z/H to bb, data  ###                               
# for score computation#                               
########################                               
### written by:      ###                               
### Bryan Caraway    ###                               
########################                               
##                    ##                                 
#                      #

##
import sys
import os
import pickle
import math
#
import deepsleepcfg as cfg
import processData  as prD 
import kinematicFit as kFit
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fun_library as lib
from fun_library import fill1e, fillne, deltaR, deltaPhi, invM, calc_mtb
np.random.seed(0)
##

def ZHbbAna(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df = kFit.retrieveData(files_, samples_, outDir_, getak8_=True)
    pt_cut = int(cfg.skim_ZHbb_dir.split('_')[-1][:3])
    for i_, key_ in enumerate(df.keys()):
        print(key_)
        fj_mask = df[key_]['ak8']['fj_lep_mask']
        #    print(j_mask)
        #
        fj_pt   = df[key_]['ak8']['pt']            [fj_mask]
        fj_phi  = df[key_]['ak8']['phi']           [fj_mask]
        fj_eta  = df[key_]['ak8']['eta']           [fj_mask]
        fj_E    = df[key_]['ak8']['E']             [fj_mask]
        sd_M    = df[key_]['ak8']['msoftdrop']     [fj_mask]
        bb_tag  = df[key_]['ak8']['btagDeepB']     [fj_mask]
        hbb_tag = df[key_]['ak8']['btagHbb']       [fj_mask]
        w_tag   = df[key_]['ak8']['deepTag_WvsQCD'][fj_mask]
        t_tag   = df[key_]['ak8']['deepTag_TvsQCD'][fj_mask]
        #print(((fj_pt >= 300) & (sd_M > 50) & (sd_M < 250)).counts)
        #plt.hist(((fj_pt >= 300) & (sd_M > 50) & (sd_M < 250)).counts)
        #plt.show()
        #exit()
        #
        fj_subjId1 = df[key_]['ak8']['subJetIdx1'] [fj_mask]
        fj_subjId2 = df[key_]['ak8']['subJetIdx2'] [fj_mask]
        #
        subj_pt    = df[key_]['ak8']['Subpt']
        subj_btag  = df[key_]['ak8']['SubbtagDeepB']
        #
        fj_sj1_pt   = subj_pt  [fj_subjId1[fj_subjId1 != -1]]
        fj_sj2_pt   = subj_pt  [fj_subjId2[fj_subjId2 != -1]]
        fj_sj1_btag = subj_btag[fj_subjId1[fj_subjId1 != -1]]
        fj_sj2_btag = subj_btag[fj_subjId2[fj_subjId2 != -1]]
        #
        lep_pt   = df[key_]['val']['Lep_pt']
        lep_eta  = df[key_]['val']['Lep_eta']
        lep_phi  = df[key_]['val']['Lep_phi']
        lep_E    = df[key_]['val']['Lep_E']
        #
        met_pt   = df[key_]['val']['MET_pt']
        met_phi  = df[key_]['val']['MET_phi']
        #
        ak8_bbcut =  ((fj_pt > pt_cut)  & (bb_tag >= 0.9))
        ak8_hbbcut = ((fj_pt >= pt_cut) & (sd_M > 50) & (sd_M < 200) & (hbb_tag >= 0.0))
        #
        tmp_ = df[key_]['df']
        b_disc = []
        pt_    = []
        phi_   = []
        eta_   = []
        E_     = []
        for key_str in tmp_.keys():
            if   ('btagDeepB_' in key_str):
                b_disc.append(key_str)
            elif ('pt_' in key_str):
                pt_.append(key_str)
            elif ('eta_' in key_str):
                eta_.append(key_str)
            elif ('phi_' in key_str):
                phi_.append(key_str)
            elif ('E_' in key_str):
                E_.append(key_str)
        b_wp = .4941
        #
        j_mask  = np.vstack(np.array(df[key_]['val']['j_lep_mask'].values))
        j_mask[j_mask == False] = True
        #
        b_disc= np.where(j_mask,tmp_[b_disc].to_numpy(),np.nan)
        b_pt  = np.where(j_mask,tmp_[pt_].to_numpy()   ,np.nan)
        b_phi = np.where(j_mask,tmp_[phi_].to_numpy()  ,np.nan)
        b_eta = np.where(j_mask,tmp_[eta_].to_numpy()  ,np.nan)
        b_E   = np.where(j_mask,tmp_[E_].to_numpy()    ,np.nan)
        #
        b_cut = (b_disc < b_wp)
        #
        j_pt  = np.where(j_mask,tmp_[pt_].to_numpy() ,np.nan)
        j_phi = np.where(j_mask,tmp_[phi_].to_numpy(),np.nan)
        j_eta = np.where(j_mask,tmp_[eta_].to_numpy(),np.nan)
        j_E   = np.where(j_mask,tmp_[E_].to_numpy()  ,np.nan)
        spher, aplan = lib.calc_SandA(
            np.append(j_pt, lep_pt.to_numpy()[:,np.newaxis], axis=1),
            np.append(j_eta, lep_eta.to_numpy()[:,np.newaxis], axis=1),
            np.append(j_phi, lep_phi.to_numpy()[:,np.newaxis], axis=1)
        )
        df[key_]['ak8']['spher'] = spher
        df[key_]['ak8']['aplan'] = aplan
        #
        b_disc[b_cut] = np.nan
        b_pt[b_cut]   = np.nan
        b_phi[b_cut]  = np.nan
        b_eta[b_cut]  = np.nan 
        b_E[b_cut]    = np.nan
        #
        j_pt [b_cut == False] = np.nan
        j_phi[b_cut == False] = np.nan
        j_eta[b_cut == False] = np.nan
        j_E  [b_cut == False] = np.nan
        #
        df[key_]['ak8']['nbbFatJets']  = bb_tag[ak8_bbcut].counts
        df[key_]['ak8']['nhbbFatJets'] = hbb_tag[ak8_hbbcut].counts
        df[key_]['ak8']['n_nonHZ_W'] = w_tag[(hbb_tag <.5) & (w_tag >= 0.8)].counts
        df[key_]['ak8']['n_nonHZ_T'] = w_tag[(hbb_tag <.5) & (t_tag >= 0.8)].counts
        #########
        hz_kinem_cut = ((fj_pt >= pt_cut) & (sd_M > 50) & (sd_M < 200) & (hbb_tag >= 0.0))
        H_hbbtag,H_pt,H_eta,H_phi,H_E,H_M,H_wtag,H_ttag,H_bbtag=lib.sortbyscore([hbb_tag    ,
                                                                                 fj_pt      ,
                                                                                 fj_eta     ,
                                                                                 fj_phi     ,
                                                                                 fj_E       ,
                                                                                 sd_M       ,
                                                                                 w_tag      , # 0.8 Med Wp
                                                                                 t_tag,
                                                                                 bb_tag    ],
                                                                                hbb_tag     ,
                                                                            hz_kinem_cut)
        hz_kinem_sj1_cut = ((fj_pt[fj_subjId1 != -1]>pt_cut) & (sd_M[fj_subjId1 != -1] > 50) & (sd_M[fj_subjId1 != -1] < 250) & (hbb_tag[fj_subjId1 != -1] >= 0.5))
        H_sj1_pt, H_sj1_btag, H_sj1_hbb = lib.sortbyscore([fj_sj1_pt  ,
                                                           fj_sj1_btag,
                                                           hbb_tag[fj_subjId1 != -1]],
                                                          hbb_tag[fj_subjId1 != -1]  ,
                                                          hz_kinem_sj1_cut)
        hz_kinem_sj2_cut = ((fj_pt[fj_subjId2 != -1]>pt_cut) & (sd_M[fj_subjId2 != -1] > 50) & (sd_M[fj_subjId2 != -1] < 250) & (hbb_tag[fj_subjId2 != -1] >= 0.5))
        H_sj2_pt, H_sj2_btag, H_sj2_hbb = lib.sortbyscore([fj_sj2_pt  ,
                                                           fj_sj2_btag,
                                                           hbb_tag[fj_subjId2 != -1]],
                                                          hbb_tag[fj_subjId2 != -1]  ,
                                                          hz_kinem_sj2_cut)
        
        # take the best tagged H/Z -> bb
        df[key_]['ak8']['H_score'] = H_hbbtag[:,0]
        df[key_]['ak8']['H_pt']    = H_pt[:,0]
        df[key_]['ak8']['H_eta']   = H_eta[:,0]
        df[key_]['ak8']['H_phi']   = H_phi[:,0]
        df[key_]['ak8']['H_E']     = H_E[:,0]
        df[key_]['ak8']['H_M']     = H_M[:,0]
        df[key_]['ak8']['H_Wscore']= H_wtag[:,0]
        df[key_]['ak8']['H_Tscore']= H_ttag[:,0]
        df[key_]['ak8']['H_bbscore']=H_bbtag[:,0]
        #
        H_sj_b12  = np.column_stack([H_sj1_btag[:,0], H_sj2_btag[:,0]])
        H_sj_pt12 = np.nan_to_num(np.column_stack([H_sj1_pt[:,0],    H_sj2_pt[:,0]]) )
        H_sjpt12_over_fjpt = (H_sj_pt12[:,0] +  H_sj_pt12[:,1])/H_pt[:,0]
        H_sjpt1_over_fjpt  = (H_sj_pt12[:,0])/H_pt[:,0]
        H_sjpt2_over_fjpt  = (H_sj_pt12[:,1])/H_pt[:,0]
        df[key_]['ak8']['n_H_sj_btag'] = np.sum(H_sj_b12 >= b_wp, axis=1)
        df[key_]['ak8']['n_H_sj']      = np.sum(H_sj_b12 >= 0, axis=1)
        df[key_]['ak8']['H_sj_bestb']  = np.nan_to_num(np.nanmax(H_sj_b12, axis=1))
        df[key_]['ak8']['H_sj_worstb'] = np.nan_to_num(np.min(H_sj_b12, axis=1))
        df[key_]['ak8']['H_sj_bbscore']= np.sum(np.nan_to_num(H_sj_b12), axis=1)
        df[key_]['ak8']['H_sjpt12_over_fjpt'] = H_sjpt12_over_fjpt
        df[key_]['ak8']['H_sjpt1_over_fjpt'] = H_sjpt1_over_fjpt
        df[key_]['ak8']['H_sjpt2_over_fjpt'] = H_sjpt2_over_fjpt
        #
        #df[key_]['ak8']['H2_score'] = H_hbbtag[:,1]
        #df[key_]['ak8']['H2_pt']    = H_pt   [:,1]
        #df[key_]['ak8']['H2_eta']   = H_eta  [:,1]
        #df[key_]['ak8']['H2_phi']   = H_phi  [:,1]
        #df[key_]['ak8']['H2_M']     = H_M    [:,1]
        #df[key_]['ak8']['H2_Wscore']= H_wtag [:,1]
        #df[key_]['ak8']['H2_bbscore']=H_bbtag[:,1]
        #########
        fjbb_pt  = fj_pt[ak8_bbcut]
        df[key_]['ak8']['fjbb_pt'] = fill1e(fjbb_pt[:,0:1] )
        df[key_]['ak8']['fjbb_score'] = fill1e(bb_tag[ak8_bbcut][:,0:1])
        fjbb_eta = fj_eta[ak8_bbcut]
        fjbb_phi = fj_phi[ak8_bbcut]
        fjbb_M   = sd_M[ak8_hbbcut]
        df[key_]['ak8']['fjbb_M'] = fill1e(fjbb_M[:,0:1])
        #
        #H_eta[H_eta.counts == 0] = np.nan
        #
        Hb_dr = deltaR(
            H_eta[:,0],H_phi[:,0],
            b_eta,b_phi)
        Hq_dr = deltaR(
            H_eta[:,0],H_phi[:,0],
            j_eta,j_phi)
        Hl_dr = deltaR(
            H_eta[:,0],H_phi[:,0],
            lep_eta,lep_phi)
        Hl_invm = invM(
            H_pt[:,0],H_eta[:,0],H_phi[:,0],
            lep_pt,lep_eta,lep_phi)
        Hl_invm_sd = lib.invM_sdM(
            H_pt[:,0],H_eta[:,0],H_phi[:,0], H_M[:,0],
            lep_pt,lep_eta,lep_phi,lep_E)
        Hl_invm_E = lib.invM_E(
            H_pt[:,0],H_eta[:,0],H_phi[:,0], H_E[:,0],
            lep_pt,lep_eta,lep_phi,lep_E)
        lb_dr = deltaR(
            lep_eta.values,lep_phi.values,
            b_eta,b_phi)
        lb_invm_E = lib.invM_E(
            lep_pt.values,lep_eta.values,lep_phi.values,lep_E.values,
            b_pt,b_eta,b_phi,b_E)
        import uproot_methods
        b_tlv = uproot_methods.TLorentzVectorArray.from_ptetaphi( b_pt.T,b_eta.T,b_phi.T,b_E.T)
        lep_tlv = uproot_methods.TLorentzVectorArray.from_ptetaphi( lep_pt,lep_eta,lep_phi,lep_E)
        bl_tlv = lep_tlv + b_tlv
        lb_mtb = lib.calc_mtb(bl_tlv.pt.T,bl_tlv.phi.T,met_pt.values,met_phi.values)            
        #
        #ind = np.argsort(np.where(Hb_dr < 0.8, Hb_dr, np.nan),axis=1)
        ind_lb = np.argsort(lb_dr,axis=1) 
        lb_invm_E_dr= np.take_along_axis(lb_invm_E,ind_lb,axis=1)
        lb_mtb_dr   = np.take_along_axis(lb_mtb,ind_lb,axis=1)
        lb_dr_dr    = np.take_along_axis(lb_dr,ind_lb,axis=1)
        max_lb_dr = np.nanmax(lb_dr_dr, axis=1)
        min_lb_dr = np.nanmin(lb_dr_dr, axis=1)
        max_lb_invm = np.nanmax(lb_invm_E_dr, axis=1)
        min_lb_invm = np.nanmin(lb_invm_E_dr, axis=1)
        df[key_]['ak8']['max_lb_dr'] = max_lb_dr
        df[key_]['ak8']['min_lb_dr'] = min_lb_dr
        df[key_]['ak8']['max_lb_invm'] = max_lb_invm
        df[key_]['ak8']['min_lb_invm'] = min_lb_invm
        #
        n_nonHbb = np.count_nonzero(Hb_dr > .8, axis=1)
        n_qnonHbb= np.count_nonzero(Hq_dr > .8, axis=1)
        n_b_Hbb  = np.count_nonzero(Hb_dr < .8, axis=1)
        n_q_Hbb  = np.count_nonzero(Hq_dr < .8, axis=1)
        #
        veto_dr = np.where(Hb_dr> 0, Hb_dr,Hq_dr)
        veto_ind= np.where(veto_dr>=.8,veto_dr,np.nan)
        rt_disc = df[key_]['valRC']['ResolvedTopCandidate_discriminator']
        rt_id   = df[key_]['valRC']['ResolvedTopCandidate_j1j2j3Idx']
        best_rt_score = []
        for idx1_, (i,j) in enumerate(zip(veto_ind,rt_id)):
            for idx2_, k in enumerate(j):                    
                if (k == '0.0.0'):
                    rt_disc[idx1_][idx2_] = np.nan
                    continue
                j1_ = int(k.split('.')[0])-1
                j2_ = int(k.split('.')[1])-1
                j3_ = int(k.split('.')[2])-1
                if ((np.isnan(i[j1_])) or (np.isnan(i[j2_])) or (np.isnan(i[j3_]))):
                    rt_disc[idx1_][idx2_] = np.nan
            #
            if (len(rt_disc[idx1_]) > 0):
                best_rt_score.append(np.nanmax(rt_disc[idx1_]))
            else:
                best_rt_score.append(np.nan)
        #
        best_rt_score = np.nan_to_num(np.array(best_rt_score))
        #print(Hb_dr)
        #print(Hj_dr)
        #exit()
        #
        ind = np.argsort(np.where(Hb_dr > 0.8, Hb_dr, np.nan),axis=1)
        b_pt_dr  = np.take_along_axis(b_pt,ind,axis=1)
        b_eta_dr = np.take_along_axis(b_eta,ind,axis=1)
        b_phi_dr = np.take_along_axis(b_phi,ind,axis=1)
        b_E_dr   = np.take_along_axis(b_E,ind,axis=1) 
        b_disc_dr= np.take_along_axis(b_disc,ind,axis=1)
        b1_tlv_dr = uproot_methods.TLorentzVectorArray.from_ptetaphi(b_pt_dr[:,0],b_eta_dr[:,0],b_phi_dr[:,0],b_E_dr[:,0])
        b2_tlv_dr = uproot_methods.TLorentzVectorArray.from_ptetaphi(b_pt_dr[:,1],b_eta_dr[:,1],b_phi_dr[:,1],b_E_dr[:,1])
        b12_pt_dr =  (b1_tlv_dr+b2_tlv_dr).pt
        Hb_invM = invM(
            H_pt[:,0],H_eta[:,0],H_phi[:,0],
            b_pt_dr,b_eta_dr,b_phi_dr)
        Hb_invM_sd = lib.invM_sdM(
            H_pt[:,0] ,H_eta[:,0] ,H_phi[:,0] ,H_M[:,0] ,
            b_pt_dr,b_eta_dr,b_phi_dr,b_E_dr)
        Hb_invM_E = lib.invM_E(
            H_pt[:,0] ,H_eta[:,0] ,H_phi[:,0] ,H_E[:,0] ,
            b_pt_dr,b_eta_dr,b_phi_dr,b_E_dr)
        mtb1 = calc_mtb(b_pt_dr[:,0],b_phi_dr[:,0],met_pt,met_phi)
        mtb2 = calc_mtb(b_pt_dr[:,1],b_phi_dr[:,1],met_pt,met_phi)
        best_Wb_invM = np.where(((mtb2 > mtb1) & (mtb2 != np.nan)), Hb_invM[:,1], Hb_invM[:,0]) 
        best_Wb_invM_sd = np.where(((mtb2 > mtb1) & (mtb2 != np.nan)), Hb_invM_sd[:,1], Hb_invM_sd[:,0]) 
        best_Wb_invM_E = np.where(((mtb2 > mtb1) & (mtb2 != np.nan)), Hb_invM_E[:,1], Hb_invM_E[:,0]) 
        #
        #print(pd.Series(Hl_dr).values.tolist())
        df[key_]['ak8']['b1_outH_score'] = b_disc_dr[:,0]
        df[key_]['ak8']['b2_outH_score'] = b_disc_dr[:,1]
        df[key_]['ak8']['b1oHZpt']       = b_pt_dr[:,0]/H_pt[:,0]
        df[key_]['ak8']['b2oHZpt']       = b_pt_dr[:,1]/H_pt[:,0]
        df[key_]['ak8']['bboHZpt']       = b12_pt_dr/H_pt[:,0]
        df[key_]['ak8']['mtb1_outH']  = mtb1
        df[key_]['ak8']['mtb2_outH']  = mtb2
        df[key_]['ak8']['best_Wb_invM']    = best_Wb_invM
        df[key_]['ak8']['best_Wb_invM_sd'] = best_Wb_invM_sd
        df[key_]['ak8']['best_Wb_invM_E']  = best_Wb_invM_E
        df[key_]['ak8']['Hb_invM1']  = Hb_invM[:,0]
        df[key_]['ak8']['Hb_invM2']  = Hb_invM[:,1]
        df[key_]['ak8']['Hb_invM1_sd']  = Hb_invM_sd[:,0]
        df[key_]['ak8']['Hb_invM2_sd']  = Hb_invM_sd[:,1]
        df[key_]['ak8']['Hb_invM1_E']   = Hb_invM_E[:,0]
        df[key_]['ak8']['Hb_invM2_E']   = Hb_invM_E[:,1]
        df[key_]['ak8']['nonHbbq1_pt'] = -np.sort(-np.where(Hq_dr > 0.8, j_pt, np.nan),axis = 1 )[:,0]
        df[key_]['ak8']['nonHbbq2_pt'] = -np.sort(-np.where(Hq_dr > 0.8, j_pt, np.nan),axis = 1 )[:,1]
        df[key_]['ak8']['nonHbbq1_dr'] = np.sort(np.where(Hq_dr > 0.8, Hq_dr, np.nan),axis = 1 )[:,0]
        df[key_]['ak8']['nonHbbq2_dr'] = np.sort(np.where(Hq_dr > 0.8, Hq_dr, np.nan),axis = 1 )[:,1]
        df[key_]['ak8']['nonHbb_b1_dr'] = np.sort(np.where(Hb_dr > 0.8, Hb_dr, np.nan),axis = 1 )[:,0]
        df[key_]['ak8']['nonHbb_b2_dr'] = np.sort(np.where(Hb_dr > 0.8, Hb_dr, np.nan),axis = 1 )[:,1]
        df[key_]['ak8']['n_nonHbb'] = n_nonHbb
        df[key_]['ak8']['n_b_Hbb'] = n_b_Hbb
        df[key_]['ak8']['n_qnonHbb'] = n_qnonHbb
        df[key_]['ak8']['n_q_Hbb'] = n_q_Hbb
        
        df[key_]['ak8']['Hl_dr']  = Hl_dr
        df[key_]['ak8']['Hl_invm']  = Hl_invm
        df[key_]['ak8']['Hl_invm_sd']  = Hl_invm_sd
        df[key_]['ak8']['Hl_invm_E']  = Hl_invm_E
        df[key_]['ak8']['best_rt_score'] = best_rt_score

        df[key_]['ak8']['lb_mtb1'] = lb_mtb_dr[:,0]
        df[key_]['ak8']['lb_invm1'] = lb_invm_E_dr[:,0]
        df[key_]['ak8']['lb_dr1'] = lb_dr_dr[:,0]
        df[key_]['ak8']['lb_mtb2'] = lb_mtb_dr[:,1]
        df[key_]['ak8']['lb_invm2'] = lb_invm_E_dr[:,1]
        df[key_]['ak8']['lb_dr2'] = lb_dr_dr[:,1]
    for key_ in df.keys():
        sample_, year_ = key_.split('_')
        with open(outDir_+'result_'+year_+'_'+sample_+'_ak8.pkl'   ,'wb') as handle:
            pickle.dump(df[key_]['ak8'], handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def plotAna(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df = kFit.retrieveData(files_, #samples_,
                           ['TTZH','TTBarLep'], 
                           outDir_, getgen_=False, getak8_=True)
    genMatched    = False
    sepGenMatched = False
    sepGen        = True
    suf = '_2017'
    df_ = {}
    if (genMatched):
        df_['TTZH_GenMatch'+suf]   = df['TTZH'+suf]
        df_['TTZH_noGenMatch'+suf] = df['TTZH'+suf]
        if ( not sepGenMatched) : 
            del df['TTZH'+suf]
    if (genMatched and sepGenMatched):                   
        df_['TTZH_genZbb'+suf]       = df['TTZH'+suf]
        df_['TTZH_genHbb'+suf]       = df['TTZH'+suf]
        df_['TTZH_genZqq'+suf]       = df['TTZH'+suf]
        del df['TTZH'+suf], df_['TTZH_GenMatch'+suf]
    #
    if (sepGen):
        df_['TTZH_Zbb'+suf]       = df['TTZH'+suf]
        df_['TTZH_Hbb'+suf]       = df['TTZH'+suf]
        #df_['TTZH_Zqq'+suf]       = df['TTZH'+suf]
        del df['TTZH'+suf] 
    if (sepGen or genMatched) : 
        df_.update(df)
        df = df_
        del df_
    print(df.keys())
    from fun_library import StackedHisto
    #StackedHisto(df, 'PrefireWeight_Down',              (-1,1), 'PrefireWeight_Down',   15)
    #StackedHisto(df, 'ISRWeight_Down',              (-2,-2), 'ISRWeight_Down',   15)
    #StackedHisto(df, 'genZHpt',              (0,450), 'genZHpt',   15)
    StackedHisto(df, 'NN',              (0,1), 'NN_output',   20)
    StackedHisto(df, 'Lep_pt',           (0,500), 'Lep_pt',    40)
    StackedHisto(df, 'H_M',     (50,200),    'HZ_M',  15)
    #StackedHisto(df, 'spher',          (0,1),  'sphericity',   20)
    #StackedHisto(df, 'aplan',          (0,.5), 'aplanarity',   20)
    #StackedHisto(df, 'nonHbb_b1_dr',    (0,5), 'nonHbb_b1_dr', 20)
    #StackedHisto(df, 'nonHbb_b2_dr',    (0,5), 'nonHbb_b2_dr', 20)
    #StackedHisto(df, 'n_q_Hbb', (0,6),      'nq_HZbb',  6) 
    #StackedHisto(df, 'H_sjpt12_over_fjpt', (0,3), 'H_sjpt12_over_fjpt', 15)
    #StackedHisto(df, 'H_sjpt1_over_fjpt', (0,3), 'H_sjpt1_over_fjpt', 15)
    #StackedHisto(df, 'H_sjpt2_over_fjpt', (0,3), 'H_sjpt2_over_fjpt', 15)
    #StackedHisto(df, 'max_lb_dr',   (0,5),     'max_lb_dr',    20)
    #StackedHisto(df, 'min_lb_dr',   (0,5),     'min_lb_dr',    20)
    #StackedHisto(df, 'max_lb_invm',   (0,750),     'max_lb_invm',    20)
    #StackedHisto(df, 'min_lb_invm',   (0,750),     'min_lb_invm',    20)
    #StackedHisto(df, 'lb_mtb1',       (0,500),    'lb_mtb1',      20)
    #StackedHisto(df, 'lb_mtb2',       (0,500),    'lb_mtb2',      20)
    #StackedHisto(df, 'lb_invm1',       (0,500),    'lb_invm1',      20)
    #StackedHisto(df, 'lb_invm2',       (0,500),    'lb_invm2',      20)
    #StackedHisto(df, 'lb_dr1',       (0,5),    'lb_dr1',      10)
    #StackedHisto(df, 'lb_dr2',       (0,5),    'lb_dr2',      10)
    #StackedHisto(df, 'n_H_sj_btag', (0,6),     'n_H_sj_btag',  6) 
    StackedHisto(df, 'nJets30',         (0,12),     'nJets30', 12)
    #StackedHisto(df, 'H_score', (.4,1),     'HZbb_score',  20)
    #StackedHisto(df, 'best_rt_score', (.5,1), 'best_rt_score', 20)
    #StackedHisto(df, 'n_qnonHbb', (0,6),     'nq_nonHZbb',  6)
    #StackedHisto(df, 'n_nonHbb', (0,6),     'nb_nonHZbb',  6)  
    #StackedHisto(df, 'Hl_dr',    (0,5),     'HZl_dr',  20)
    #StackedHisto(df, 'Hl_invm',  (0,700),   'HZl_invm',  50)
    #StackedHisto(df, 'Hl_invm_sd',  (0,700),   'HZl_invm_sd',  50)
    #StackedHisto(df, 'Hl_invm_E',  (0,700),   'HZl_invm_E',  50)
    #StackedHisto(df, 'n_H_sj', (0,4),     'n_H_sj',  4) 
    #StackedHisto(df, 'n_b_Hbb', (0,6),      'nb_HZbb',  6) 
    #StackedHisto(df, 'H_sj_bestb', (0,1), 'H_sj_bestb', 20)
    #StackedHisto(df, 'H_sj_worstb', (0,1), 'H_sj_worstb', 20)
    #StackedHisto(df, 'H_sj_bbscore', (0,2), 'H_sj_bbscore', 20)
    #StackedHisto(df, 'nMergedTops', (0,5),'nMergedTops',5)
    #StackedHisto(df, 'n_nonHZ_W', (0,4),       'n_nonHZ_W', 4)
    #StackedHisto(df, 'n_nonHZ_T', (0,4),       'n_nonHZ_T', 4)
    StackedHisto(df, 'H_eta',     (-3.2,3.2),    'HZ_eta',  50)
    #StackedHisto(df, 'H_bbscore', (0,1), 'H_bbscore', 20)
    #StackedHisto(df, 'b1_outH_score', (0,1), 'b1_outH_score', 20)
    #StackedHisto(df, 'b2_outH_score', (0,1), 'b2_outH_score', 20)
    #StackedHisto(df, 'b1oHZpt',       (0,5), 'b1oHZpt',       20)
    #StackedHisto(df, 'b2oHZpt',       (0,5), 'b2oHZpt',       20)
    #StackedHisto(df, 'bboHZpt',       (0,5), 'bboHZpt',       20)
    #StackedHisto(df, 'best_Wb_invM',    (0,1000), 'best_Wb_invM',    50)
    #StackedHisto(df, 'Hb_invM1',        (0,1000), 'Hb_invM1',        50)
    #StackedHisto(df, 'Hb_invM2',        (0,1000), 'Hb_invM2',        50)
    #StackedHisto(df, 'best_Wb_invM_sd', (0,1000), 'best_Wb_invM_sd', 50)
    #StackedHisto(df, 'Hb_invM1_sd',     (0,1000), 'Hb_invM1_sd',     50)
    #StackedHisto(df, 'Hb_invM2_sd',     (0,1000), 'Hb_invM2_sd',     50)
    #StackedHisto(df, 'best_Wb_invM_E',  (0,1000), 'best_Wb_invM_E',  50)
    #StackedHisto(df, 'Hb_invM1_E',      (0,1000), 'Hb_invM1_E',      50)
    #StackedHisto(df, 'Hb_invM2_E',      (0,1000), 'Hb_invM2_E',      50)
    #StackedHisto(df, 'H2_M',     (0,300),    'HZ2_M',  40)
    StackedHisto(df, 'H_pt',     (200,600), 'HZ_pt',  20)
    #StackedHisto(df, 'H2_pt',     (200,600), 'HZ2_pt',  20)
    #StackedHisto(df, 'H2_score', (-1,1),     'HZbb2_score',  20)
    StackedHisto(df, 'H_Wscore', (0,1),     'H_Wscore',  20)
    StackedHisto(df, 'H_Tscore', (0,1),     'H_Tscore',  20)
   
    StackedHisto(df, 'MET_pt',  (20,500),   'MET',         40)
    
    StackedHisto(df, 'mtb1_outH',  (0,500), 'mtb1_outH',  40)
    StackedHisto(df, 'mtb2_outH',  (0,500), 'mtb2_outH',  40)
    StackedHisto(df, 'nhbbFatJets', (0,6),      'nhbbFatJets',  6) 
    
    #StackedHisto(df, 'H2_Wscore', (0,1),     'H2_Wscore',  20)
    StackedHisto(df, 'nFatJets',      (0,5),     'nFatJets', 5)
    StackedHisto(df, 'nJets',         (0,12),     'nJets', 12)
    #StackedHisto(df, 'nonHbbj1_pt',     (0,600), 'nonHbbq1_pt',  60)
    #StackedHisto(df, 'nonHbbj2_pt',     (0,600), 'nonHbbq2_pt',  60)
    #StackedHisto(df, 'nonHbbj1_dr',     (0,5), 'nonHbbq1_dr',  20)
    #StackedHisto(df, 'nonHbbj2_dr',     (0,5), 'nonHbbq2_dr',  20)
    
    StackedHisto(df, 'nResolvedTops', (0,5),'nResolvedTops',5)

def Bkg_Est(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    from fun_library import calc_Kappa
    df = kFit.retrieveData(files_, samples_, outDir_, getgen_=False, getak8_=True)
    sig = 'TTZH'
    bkg = 'TTBarLep'
    suf = '_2017'
    #
    steps = 0.05
    n_bins = 13
    #
    bins = np.arange(0.0,.95,steps)
    m_bins = np.arange(50.0,200.0+(200.0-50.0)/n_bins,(200.0-50.0)/n_bins)
    kap = kap_e = c_eff = c_eff_e = np.zeros(len(bins))
    for i_ in range(len(bins)):
        #kap[i_], kap_e[i_], c_eff[i_], c_eff_e[i_] = calc_Kappa(df,[bins[i_],bins[i_]+steps])
        #bkg, bkg_err, eff, eff_err = calc_Kappa(df,[bins[i_],bins[i_]+steps]) 
        #if i_ == 0:
        #    plt.errorbar(x=[60,85,122.5,172.5], y=eff, xerr=[10,15,22.5,27.5], yerr=eff_err,
        #                 fmt='o',barsabove=True,capsize=5)
        #    plt.title('Sig/Bkg for NN score > .95')
        #    plt.show()
        #    plt.clf()
        #plt.errorbar(x=[60,85,122.5,172.5], y=bkg, xerr=[10,15,22.5,27.5], yerr=bkg_err,
        #             fmt='o',barsabove=True,capsize=5)
        #plt.title('Bkg ratio for NN score range: {:.2f}-->{:.2f}'.format(bins[i_],bins[i_]+steps))
        #plt.show()
        #plt.clf()
        def applyCuts(df_,key_):
            base_cuts =(
                (df_[key_]['ak8']['n_nonHbb'] >= 2)    &
                (df_[key_]['ak8']['nhbbFatJets'] > 0)  &
                (df_[key_]['ak8']['H_M']         > 50) &
                (df_[key_]['ak8']['H_M']         < 200))
            if (sig in key_): # look at the sig/bkg ratios only in the signal region
                base_cuts = base_cuts & (df_[key_]['val']['NN'] >= .95)
                _genm   = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == True)
                _nogenm = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == False)
                _genZ   = base_cuts & (df_[key_]['val']['matchedGen_Zbb'] == True)
                _genH   = base_cuts & (df_[key_]['val']['matchedGen_Hbb'] == True)
                sig_genm, sig_genmW     = df_[key_]['ak8']['H_M'][_genm],   (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genm  ]
                sig_nogenm, sig_nogenmW = df_[key_]['ak8']['H_M'][_nogenm], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_nogenm]
                sig_genZ, sig_genZW     = df_[key_]['ak8']['H_M'][_genZ],   (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genZ  ]
                sig_genH, sig_genHW     = df_[key_]['ak8']['H_M'][_genH],   (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genH  ]

                return [sig_genm, sig_nogenm, sig_genZ, sig_genH,
                        sig_genmW, sig_nogenmW, sig_genZW, sig_genHW]
            else:
                _cr = base_cuts & (df_[key_]['val']['NN'] >= bins[i_]) & (df_[key_]['val']['NN'] < bins[i_]+steps)
                bkg_cr, bkg_cr_gw = df_[key_]['ak8']['H_M'][_cr], np.sign(df_[key_]['val']['genWeight'][_cr])
                _sr = base_cuts & (df_[key_]['val']['NN'] >= .95)
                bkg_sr, bkg_srW, bkg_sr_gw = df_[key_]['ak8']['H_M'][_sr], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_sr  ], np.sign(df_[key_]['val']['genWeight'][_sr])  
                
                return bkg_cr, bkg_sr, bkg_srW, bkg_cr_gw, bkg_sr_gw
        #
        def shapePlot(y_,y_err_,label_,ylabel_):
            bin_c = (m_bins[1:] + m_bins[:-1])/2
            bin_w = m_bins[1:] - m_bins[:-1]
            plt.errorbar(x=bin_c, y=y_, xerr=bin_w/2, yerr=y_err_,
                         fmt='o',barsabove=True,capsize=5, label=label_) 
            #plt.title(title_)
            plt.xlabel('M_SD')
            plt.ylabel(ylabel_)
            plt.ylim(0,None)
            #plt.show()
            #plt.clf()
        def normed(x_,val_):
            bin_w = m_bins[1:] - m_bins[:-1]
            n_events = float(sum(val_))
            return x_/n_events/bin_w
        def calc_ratioE(a_,b_,a_w=None, b_w=None, norm=False):
            a_hist,_ = np.histogram(a_,m_bins,weights=a_w)
            b_hist,_ = np.histogram(b_,m_bins,weights=b_w)
            c_= a_hist/b_hist
            if ((abs(a_w) != 1).any() or (abs(b_w) != 1).any()):
                a_sumw2,_ = np.histogram(a_,m_bins,weights=np.power(a_w,2))
                b_sumw2,_ = np.histogram(b_,m_bins,weights=np.power(b_w,2))
                a_err, b_err = np.sqrt(a_sumw2), np.sqrt(b_sumw2)
            else:
                a_err, b_err = np.sqrt(a_hist), np.sqrt(b_hist)
            if norm:
                a_err, b_err   = normed(a_err,a_hist), normed(b_err,b_hist)
                a_hist, b_hist = normed(a_hist,a_hist), normed(b_hist,b_hist)
            c_err_ = abs(c_)*np.sqrt(np.power(a_err/a_hist,2)+np.power(b_err/b_hist,2))
            return c_,c_err_
        def histvals(vals_,val_w=None,norm=False):
            hist_vals,_ =  np.histogram(vals_,m_bins,weights=val_w)
            val_sumw2,_ = np.histogram(vals_,m_bins,weights=np.power(val_w,2))
            hist_vals_err = np.sqrt(val_sumw2)
            if norm:
                return normed(hist_vals,hist_vals), normed(hist_vals_err,hist_vals)
            return hist_vals, hist_vals_err
        #
        bkg_cr, bkg_sr, bkg_srW, bkg_cr_gw, bkg_sr_gw = applyCuts(df,bkg+suf)
        bkg_df = pd.DataFrame(columns=['bkg_sr','bkg_srW'])
        for key_ in df.keys():
            if sig in key_: continue
            _cr, _sr, _srW, cr_gw, sr_gw = applyCuts(df,key_)
            if _sr.size == 0 : continue            
            bkg_df = bkg_df.append(pd.DataFrame({'bkg_sr':_sr,'bkg_srW':_srW}),  ignore_index = True) 
        #bkg_df = pd.DataFrame.from_records([applyCuts(df,x) if sig not in x else np.array([[],[],[],[],[]]) for x in df.keys()], columns=['cr','sr','srW','cr_gw','sr_gw'])
        if (i_ >= 0):
            sig_genm, sig_nogenm, sig_genZ, sig_genH, sig_genmW, sig_nogenmW, sig_genZW, sig_genHW = applyCuts(df,sig+suf)
            #shapePlot(*calc_ratioE(sig_genm,bkg_sr,  a_w=sig_genmW,    b_w=bkg_srW, norm=True),    'Sig/Bkg SR(genM)',   'Sig/Bkg SR(genM)')
            #shapePlot(*calc_ratioE(sig_nogenm,bkg_sr,a_w=sig_nogenmW,  b_w=bkg_srW, norm=True),'Sig/Bkg SR(nogenM)', 'Sig/Bkg SR(nogenM)')
            #shapePlot(*calc_ratioE(sig_genZ,bkg_sr,  a_w=sig_genZW,    b_w=bkg_srW, norm=True),    'Sig/Bkg SR(genZ)',   'Sig/Bkg SR(genZ)')
            #shapePlot(*calc_ratioE(sig_genH,bkg_sr,  a_w=sig_genHW,    b_w=bkg_srW, norm=True),    'Sig/Bkg SR(genH)',   'Sig/Bkg SR(genH)')
            # For Ken: Shape in Sig Region
            # 1: Shape of ttH, ttZ, tt
            shapePlot(*histvals(bkg_df['bkg_sr'],  val_w=bkg_df['bkg_srW'],   norm=False),    'bkg',   '')    
            shapePlot(*histvals(sig_genH,          val_w=sig_genHW,        norm=False),    'ttH',   '')    
            shapePlot(*histvals(sig_genZ,          val_w=sig_genZW,        norm=False),    'ttZ',   '')    

            plt.legend()
            plt.show()
            plt.clf()
            # 2: Shape of genM vs noGenM
            shapePlot(*histvals(sig_genm,   val_w=sig_genmW,   norm=False),    'norm(GenMatch)',     '')    
            shapePlot(*histvals(sig_nogenm, val_w=sig_nogenmW, norm=False),    'norm(noGenMatch)',   '')    
            plt.legend()
            plt.show()
            plt.clf()
        #
        shapePlot(*calc_ratioE(bkg_sr,bkg_cr, a_w=bkg_sr_gw, b_w=bkg_cr_gw, norm=True), 'BKG_SR/BKG_CR NN score range {:.2f}-->{:.2f}'.format(bins[i_],bins[i_]+steps), 'BKG_SR/BKG_CR')


        
        

def GenAna_ttbar(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df_ = kFit.retrieveData(files_, ['TTBarLep'], outDir_, getgen_=True, getak8_=True)
    for key_ in df_.keys():
        df = df_[key_]
        w  = df['val']['weight'].values * np.sign(df['val']['genWeight'].values) * (137/41.9)
        #
        gen_df = df['gen']
        fat_df = df['ak8']
        met    = df['val']['MET_pt']
        base_cuts = (
            (fat_df['n_nonHbb']   >=  2) &  
            #(fat_df['n_b_Hbb']    >=  1) &
            (fat_df['nhbbFatJets']>   0) & 
            (fat_df['H_M']        >  50) & 
            (fat_df['H_M']        < 200) &
            #(fat_df['best_Wb_invM']> 200) &
            #(fat_df['H_Wscore'] < .90) &
            (met                  >=  0)
        )
        #w = w[base_cuts]
        #
        ZH_eta   = fat_df['H_eta']           #[base_cuts]
        ZH_phi   = fat_df['H_phi']           #[base_cuts]
        ZH_M     = fat_df['H_M']             #[base_cuts]
        ZH_score = fat_df['H_score']         #[base_cuts] 
        ZH_Wscore= fat_df['H_Wscore']        #[base_cuts]
        best_Wb_invM = fat_df['best_Wb_invM']#[base_cuts]
        #
        #fig, ax = plt.subplots()
        #ax.hist2d(x=ZH_Wscore, y=best_Wb_invM, 
        #          range= ((0,1),(0,300)),
        #          cmin = 0.01,
        #          bins=50, weights=w)
        ##plt.show()
        #plt.clf()
        #
        gen_ids = gen_df['GenPart_pdgId']           #[base_cuts]
        gen_mom = gen_df['GenPart_genPartIdxMother']#[base_cuts]
        gen_pt  = gen_df['GenPart_pt']              #[base_cuts]
        gen_eta = gen_df['GenPart_eta']             #[base_cuts]
        gen_phi = gen_df['GenPart_phi']             #[base_cuts]
        gen_E   = gen_df['GenPart_E']               #[base_cuts]
        #
        
        get_bq_fromtop = (
            ((abs(gen_ids) == 5) & (abs(gen_ids[gen_mom]) == 6)) | # if bottom and mom is top
            ((abs(gen_ids) < 5) & (abs(gen_ids[gen_mom]) == 24) & (abs(gen_ids[gen_mom[gen_mom]]) == 6))   # get quarks from W whos parrent is a top
        )
        base = ((get_bq_fromtop) & (get_bq_fromtop.sum() ==4))
        case_1_base_pt = ((base) & 
                          (gen_ids[gen_mom] > 0) & (gen_ids[gen_mom] != 21))
        case_1_base_mt = ((base) & 
                          (gen_ids[gen_mom] < 0) & (gen_ids[gen_mom] != -21))
        case_1_base    = (
            ((case_1_base_pt) & (case_1_base_pt.sum() == 3)) |
            ((case_1_base_mt) & (case_1_base_mt.sum() == 3)) 
        )
        case_4_base = ((base & (abs(gen_ids) == 5)) & ((base & (abs(gen_ids) == 5)).sum() == 2))
        case_5_base_npt = (base & ((gen_ids[gen_mom] == -6) | (gen_ids[gen_mom] ==  24)))
        case_5_base_nmt = (base & ((gen_ids[gen_mom] ==  6) | (gen_ids[gen_mom] == -24)))
        case_5_base = (
            ((case_5_base_npt) & (case_5_base_npt.sum() == 3)) |
            ((case_5_base_nmt) & (case_5_base_nmt.sum() == 3)) 
        )   
        # 
        case_123_dr = deltaR(ZH_eta,ZH_phi,gen_eta[case_1_base],gen_phi[case_1_base])
        case_4_dr   = deltaR(ZH_eta,ZH_phi,gen_eta[case_4_base],gen_phi[case_4_base])
        case_56_dr   = deltaR(ZH_eta,ZH_phi,gen_eta[case_5_base],gen_phi[case_5_base])
        #
        case_1  = (((case_123_dr < 0.8) & 
                    (((abs(gen_ids[case_1_base][case_123_dr<0.8]) < 6).sum() == 2) & ((abs(gen_ids[case_1_base][case_123_dr<0.8]) < 5).sum() == 2))
                ).sum() == 2)
        case_2  = (((case_123_dr < 0.8) & 
                    (((abs(gen_ids[case_1_base][case_123_dr<0.8]) < 6).sum() == 2) & ((abs(gen_ids[case_1_base][case_123_dr<0.8]) < 5).sum() == 1))
                ).sum() == 2)
        case_3  = (((case_123_dr < 0.8) & 
                    (((abs(gen_ids[case_1_base][case_123_dr<0.8]) < 6).sum() == 3) & ((abs(gen_ids[case_1_base][case_123_dr<0.8]) < 5).sum() == 2))
                ).sum() == 3)
        case_4  = (((case_4_dr < 0.8) & ((abs(gen_ids[case_4_base][case_4_dr<0.8]) == 5).sum() == 2)).sum() == 2)
        case_5  = (((case_56_dr < 0.8) &
                    (((abs(gen_ids[case_5_base][case_56_dr<0.8]) < 6).sum() == 2) & ((abs(gen_ids[case_5_base][case_56_dr<0.8]) < 5).sum() == 1))
                ).sum() == 2)
        case_6  = (((case_56_dr < 0.8) &
                    (((abs(gen_ids[case_5_base][case_56_dr<0.8]) < 6).sum() == 3) & ((abs(gen_ids[case_5_base][case_56_dr<0.8]) < 5).sum() == 2))
                ).sum() == 3)
        case_7  = ((case_1 == False) & (case_2 == False) & (case_3 == False) & (case_4 == False) & (case_5 == False) & (case_6 == False))
        #
        df_[key_]['val']['case_1'] = case_1
        df_[key_]['val']['case_2'] = case_2
        df_[key_]['val']['case_3'] = case_3
        df_[key_]['val']['case_4'] = case_4
        df_[key_]['val']['case_5'] = case_5
        df_[key_]['val']['case_6'] = case_6
        df_[key_]['val']['case_7'] = case_7
        year   = key_[key_.index('201'):key_.index('201')+4]
        sample = key_.split('_201')[0]
        df_[key_]['val'].to_pickle(outDir_+'result_'+year+'_'+sample+'_val.pkl')
        exit()
        #
        case_1 = case_1 & base_cuts
        case_2 = case_2 & base_cuts
        case_3 = case_3 & base_cuts
        case_4 = case_4 & base_cuts
        case_5 = case_5 & base_cuts
        case_6 = case_6 & base_cuts
        case_7 = case_7 & base_cuts
        #
        case_7_dr = deltaR(ZH_eta[case_7],ZH_phi[case_7],gen_eta[case_7],gen_phi[case_7])
        ind_ = 1
        T = len(gen_ids[case_7])
        G = (((case_7_dr < .8) & (gen_ids[case_7] ==21)).sum() > 0)
        L = (((case_7_dr < .8) & (abs(gen_ids[case_7]) >= 11) & (abs(gen_ids[case_7]) <= 16)).sum() > 0)
        B = (((case_7_dr < .8) & (abs(gen_ids[case_7]) == 5)).sum() > 0)
        BfromG = (((case_7_dr < .8) & (abs(gen_ids[case_7]) == 5) & (gen_ids[gen_mom][case_7] == 22)).sum() > 0)
        Q = (((case_7_dr < .8) & (abs(gen_ids[case_7]) < 5)).sum() > 0)
        QfromG = (((case_7_dr < .8) & (abs(gen_ids[case_7]) < 5) & (gen_ids[gen_mom][case_7] == 22)).sum() > 0)
        #
        BnG = (B & ~G)
        BnL = (B & ~L)
        BnLnG = (B & ~L & ~G)
        GnB = (G & ~B)
        GnL = (G & ~L)
        GnLnB = (G & ~L & ~B)
        LnB = (L & ~B)
        LnG = (L & ~G)
        LnGnB = (L & ~G & ~B)
        nBnGnL = (~B & ~L & ~G)
        #
        LBnG   = (L & B & ~G)
        GBnL   = (G & B & ~L)
        LGnB   = (L & G & ~B) 
        #
        GnQ    = (G & ~Q)
        LnQ    = (L & ~Q)
        BnQ    = (B & ~Q)
        print(T)
        print(G.sum(),     'G')
        print(L.sum(),     'L')
        print(B.sum(),     'B')
        print(Q.sum(),     'Q')
        print(BnG.sum(),   'BnG')
        print(BnL.sum(),   'BnL')
        print(GnL.sum(),   'GnL')
        print(GnB.sum(),   'GnB')
        print(LnB.sum(),   'LnB')
        print(LnG.sum(),   'LnG')
        print(BnLnG.sum(), 'BnLnG')
        print(GnLnB.sum(), 'GnLnB')
        print(LnGnB.sum(), 'LnGnB')
        print(nBnGnL.sum(),'nBnGnL')
        print(LBnG.sum(),  'LBnG')
        print(GBnL.sum(),  'GBnL')
        print(LGnB.sum(),  'LGnB')
        print(GnQ.sum(),   'GnQ')
        print(LnQ.sum(),   'LnQ')
        print(BnQ.sum(),   'BnQ')
        print(BfromG.sum(),'BfromG')
        print(QfromG.sum(),'QfromG')
        #print(gen_ids[case_7][(case_7_dr < .8)][B & Q])
        #print(gen_ids[gen_mom][case_7][(case_7_dr < .8)][B & Q])
        #print(case_7_dr[(case_7_dr < .8)][B & Q])
        #print(gen_ids[case_7][(case_7_dr < .8)][((case_7_dr < .8) & (abs(gen_ids[case_7]) >= 11) & (abs(gen_ids[case_7]) <= 16)).sum() > 0])
        #print(gen_ids[case_7][((case_7_dr < .8) & (gen_ids[case_7] ==21)).sum() > 0].flatten())
        plt.hist(gen_ids[case_7][(case_7_dr < .8)][G].flatten(), bins=100, range=(-25,25))
        plt.title('gluon in H/Z ({0:2.2f}%)'.format(G.sum()/T * 100))
        plt.show()
        plt.clf()
        plt.hist(gen_ids[case_7][(case_7_dr < .8)][L].flatten(), bins=100,range=(-25,25))
        plt.title('lepton in H/Z ({0:2.2f}%)'.format(L.sum()/T * 100))
        plt.show()
        plt.clf()
        plt.hist(gen_ids[case_7][(case_7_dr < .8)][B].flatten(), bins= 100, range=(-25,25))
        plt.title('b in H/Z ({0:2.2f}%)'.format(B.sum()/T * 100))
        plt.show()
        plt.clf()
        plt.hist(gen_ids[case_7][(case_7_dr < .8)][nBnGnL].flatten(), bins= 100, range=(-25,25))
        plt.title('no b,l,g in H/Z ({0:2.2f}%)'.format(nBnGnL.sum()/T * 100))
        plt.show()
        plt.clf()
        #print(gen_ids[case_7][(((case_7_dr < .8) & (gen_ids[case_7] ==21)).sum() > 0)][ind_])
        #print(*gen_ids[gen_mom][case_7][(((case_7_dr < .8) & (gen_ids[case_7] ==21)).sum() > 0)][ind_])
        #print(*case_7_dr[(((case_7_dr < .8) & (gen_ids[case_7] ==21)).sum() > 0)][ind_])
        #exit()
        #print(((case_2 == case_4) & (case_2 == True)).sum()) ########
        #print(((case_3 == case_4) & (case_3 == True)).sum()) ########
        #print(gen_ids[((case_2 == case_4) & (case_2 == True))]) ########
        #print(*gen_ids[((case_3 == case_4) & (case_3 == True))][0]) ########
        #print(*gen_ids[gen_mom][((case_3 == case_4) & (case_3 == True))][0]) ########
        cases = [case_1,case_2,case_3,case_4,case_5,case_6,case_7]
        des_dict = {
            "Case_1" : 'Case_1 (qq from W within fatjet)',
            "Case_2" : 'Case_2 (b+q from W, same top within fatjet)',
            "Case_3" : 'Case_3 (b+qq from W, same top within fatjet)',
            "Case_4" : 'Case_4 (bb from ttbar within fatjet)',
            "Case_5" : 'Case_5 (b from one top, and q from other top within fatjet)',
            "Case_6" : 'Case_6 (b from one top, and qq from other top within fatjet)',
            "case_7" : 'case_7 (Else)'
        }
        for case, des_key in zip(cases,des_dict):
            n,b,_ =plt.hist(ZH_M[case],weights=w[case],bins = 50)    
            plt.title(des_dict[des_key]+' ({0:3.1f})'.format(sum(n[:])))
            plt.xlabel('reco_HZ_M (softdrop)')
            plt.show()
            
        for case, des_key in zip(cases,des_dict):
            n,b,_ =plt.hist(ZH_score[case],weights=w[case],bins = 50)    
            plt.title(des_dict[des_key]+' ({0:3.1f})'.format(sum(n[:])))
            plt.xlabel('reco_HZ_score')
            plt.show()
        for case, des_key in zip(cases,des_dict):
            n,b,_ =plt.hist(ZH_Wscore[case],weights=w[case],bins = 50)    
            plt.title(des_dict[des_key]+' ({0:3.1f})'.format(sum(n[:])))
            plt.xlabel('reco_HZ_Wscore')
            plt.show()
        for case, des_key in zip(cases,des_dict):
            n,b,_ =plt.hist(best_Wb_invM[case],weights=w[case],bins = 50)    
            plt.title(des_dict[des_key]+' ({0:3.1f})'.format(sum(n[:])))
            plt.xlabel('reco_Wb_invM')
            plt.show()

def GenAna(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df = kFit.retrieveData(files_, ['TTZH'], outDir_, getgen_=True, getak8_=True)
    df = df['TTZH_2017']
    w  = df['val']['weight'].values * np.sign(df['val']['genWeight'].values) * (137/41.9)
    pt_cut = int(cfg.skim_ZHbb_dir.split('_')[-1][:3])
    #
    gen_df = df['gen']
    fat_df = df['ak8']
    val_df = df['val']
    #
    zh_ids = gen_df['GenPart_pdgId']
    zh_mom = gen_df['GenPart_genPartIdxMother']
    zh_pt  = gen_df['GenPart_pt']
    zh_eta = gen_df['GenPart_eta']
    zh_phi = gen_df['GenPart_phi']
    zh_E   = gen_df['GenPart_E']
    # higgs pdgid = 25
    # Z     pdgid = 23
    isHbb = (zh_ids == 25)
    #
    isbb_fromZ = (((zh_ids == -5) | (zh_ids == 5)) & 
                 (zh_ids[zh_mom] == 23))
    isZbb  = ((zh_ids == 23) & (isbb_fromZ.sum() > 0))
    isZHbb = (isHbb | isZbb)
    #
    zh_pt = fill1e(zh_pt[isZHbb]).flatten()
    zh_eta = fill1e(zh_eta[isZHbb]).flatten()
    zh_phi = fill1e(zh_phi[isZHbb]).flatten()
    zh_E   = fill1e(zh_E[isZHbb]).flatten()    
    #
    fj_pt   = fillne(fat_df['pt'])
    fj_phi  = fillne(fat_df['phi'])
    fj_eta  = fillne(fat_df['eta'])
    fj_E    = fillne(fat_df['E'])
    sd_M    = fillne(fat_df['msoftdrop'])
    fj_w_tag= fillne(fat_df['deepTag_WvsQCD'])
    fj_b_tag= fillne(fat_df['btagDeepB'])
    hbb_tag = fillne(fat_df['btagHbb'])
    best_wb = fat_df['best_Wb_invM']
    #
    gen_dr = deltaR(zh_eta,zh_phi,fj_eta,fj_phi)
    gen_dr_match = np.nanmin(gen_dr,axis=1)
    gen_match    = ((gen_dr_match < 0.8) & (zh_pt >= pt_cut)) #& (zh_eta <= 2.6) & (zh_eta >= -2.6))
    #
    def matchkinem(kinem_):
        ind_=np.argsort(gen_dr,axis=1)
        return np.take_along_axis(kinem_,ind_,axis=1)[:,0]
    reco_zh_pt    = matchkinem(fj_pt)[gen_match]
    reco_zh_M     = matchkinem(sd_M)[gen_match]
    reco_zh_score = matchkinem(hbb_tag)[gen_match]
    reco_zh_Wscore= matchkinem(fj_w_tag)[gen_match]
    reco_zh_bscore= matchkinem(fj_b_tag)[gen_match]
    reco_zh_wb_invM  = best_wb[gen_match]
    reco_zh_w     = w[gen_match]
    #isZ = (isZ.counts == 1)[gen_match]
    gen_matchZ = (gen_match & (isZbb.sum() == 1))
    reco_z_pt    = matchkinem(fj_pt)[gen_matchZ]
    reco_z_M     = matchkinem(sd_M)[gen_matchZ]
    reco_z_score = matchkinem(hbb_tag)[gen_matchZ]
    reco_z_Wscore= matchkinem(fj_w_tag)[gen_matchZ]
    reco_z_bscore= matchkinem(fj_b_tag)[gen_matchZ]
    reco_z_wb_invM= best_wb[gen_matchZ]
    reco_z_w     = w[gen_matchZ]
    #
    gen_matchH = (gen_match & (isHbb.sum() == 1))
    reco_h_pt    = matchkinem(fj_pt)[gen_matchH]
    reco_h_M     = matchkinem(sd_M)[gen_matchH]
    reco_h_score = matchkinem(hbb_tag)[gen_matchH]
    reco_h_Wscore= matchkinem(fj_w_tag)[gen_matchH]
    reco_h_bscore= matchkinem(fj_b_tag)[gen_matchH]
    reco_h_wb_invM  = best_wb[gen_matchH]
    reco_h_w     = w[gen_matchH]
    #
    # look at the ratio of events at higher eta
    # look at the efficiency for Z/H tagging at eta > 2.4
    print('\n H/Z >=',pt_cut)
    print('# of H eta total: ', len(zh_eta[(isHbb.sum() == 1) & (zh_pt >= pt_cut)]))
    print('# of H |eta| >= 2.4: ', len(zh_eta[(isHbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4)]))
    print('# of H |eta| >= 2.4, |eta| <= 2.5: ', len(zh_eta[(isHbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4) & (abs(zh_eta) <= 2.5)]))
    print('# of H |eta| >= 2.4, |eta| <= 2.5, dr matched to fatjet: ', len(zh_eta[(isHbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4) & (abs(zh_eta) <= 2.5) & (gen_match)]))
    print('# of H |eta| >= 2.4, |eta| <= 2.5, dr matched to fatjet, NN >= 0: ', len(zh_eta[(isHbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4) & (abs(zh_eta) <= 2.5) & (gen_match) & (val_df['NN'] >= 0)]))
    print('\n')
    print('# of Z eta total: ', len(zh_eta[(isZbb.sum() == 1) & (zh_pt >= pt_cut)]))
    print('# of Z |eta| >= 2.4: ', len(zh_eta[(isZbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4)]))
    print('# of Z |eta| >= 2.4, |eta| <= 2.5: ', len(zh_eta[(isZbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4) & (abs(zh_eta) <= 2.5)]))
    print('# of Z |eta| >= 2.4, |eta| <= 2.5, dr matched to fatjet: ', len(zh_eta[(isZbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4) & (abs(zh_eta) <= 2.5) & (gen_match)]))
    print('# of Z |eta| >= 2.4, |eta| <= 2.5, dr matched to fatjet, NN >= 0: ', len(zh_eta[(isZbb.sum() == 1) & (zh_pt >= pt_cut) & (abs(zh_eta) >= 2.4) & (abs(zh_eta) <= 2.5) & (gen_match) & (val_df['NN'] >= 0)]))
    exit()
    #
    from matplotlib.ticker import AutoMinorLocator
    def simpleplot(x_,l_,w_):
        for x,l in zip(x_,l_):
            fig, ax = plt.subplots()
            ax.hist(x, bins=50, weights=w_)
            plt.xlabel(l)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.grid(True)
            #plt.show()
        fig, ax = plt.subplots()
        ax.hist2d(x=x_[2], y=x_[1], 
                  range= ((-1,1),(0,300)),
                  cmin = 0.01,
                  bins=50, weights=w_)
        plt.xlabel(l_[2])
        plt.ylabel(l_[1])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True)
        #plt.show()
        #
    simpleplot([reco_zh_pt,reco_zh_M,reco_zh_score,reco_zh_Wscore,reco_zh_bscore,reco_zh_wb_invM],['zh_pt','zh_M','zh_score','reco_zh_Wscore','reco_zh_bscore','reco_zh_wb_invM'],reco_zh_w)
    simpleplot([reco_z_pt, reco_z_M, reco_z_score,reco_z_Wscore,reco_z_bscore,reco_z_wb_invM],  ['z_pt', 'z_M', 'z_score','reco_z_Wscore','reco_z_bscore','reco_z_wb_invM'],reco_z_w)
    simpleplot([reco_h_pt, reco_h_M, reco_h_score,reco_h_Wscore,reco_h_bscore,reco_h_wb_invM],  ['h_pt', 'h_M', 'h_score','reco_h_Wscore','reco_h_bscore','reco_h_wb_invM'],reco_h_w)
    #
    import matplotlib.backends.backend_pdf as matpdf
    pdf = matpdf.PdfPages('money_pdf/money_genZH.pdf')
    for fig_ in range(1, plt.gcf().number+1):
        pdf.savefig( fig_ )
    pdf.close()


def matchLep(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df_ = kFit.retrieveData(files_, ['TTZH','TTBarLep'], outDir_, getak8_=True, getgen_=True)
    for key_ in df_.keys():
        year = key_.split('_')[1]
        sample = key_.split('_')[0]
        #
        df = df_[key_]
        val_df = df['val']
        gen_df = df['gen']
        fj_df  = df['ak8']
        weight = df['val']['weight']*np.sign(df['val']['genWeight'])*(137/41.9)*\
                 df_[key_]['val']['BTagWeight']*df_[key_]['val']['puWeight']*df_[key_]['val']['PrefireWeight']
        #
        gen_ids = gen_df['GenPart_pdgId']
        gen_mom = gen_df['GenPart_genPartIdxMother']
        gen_pt  = gen_df['GenPart_pt'] 
        gen_eta = gen_df['GenPart_eta']
        gen_phi = gen_df['GenPart_phi']
        #
        lep_pt   = val_df['Lep_pt'].values
        lep_eta  = val_df['Lep_eta'].values
        lep_phi  = val_df['Lep_phi'].values
        #
        fj_pt    = fj_df['pt']
        fj_eta   = fj_df['eta']
        fj_phi   = fj_df['phi']
        fj_Hscore= fj_df['btagHbb']
        #
        lep_cut = (((abs(gen_ids) == 11) | (abs(gen_ids) == 13) | (abs(gen_ids) == 15)) & ((abs(gen_ids[gen_mom]) == 24) & (abs(gen_ids[gen_mom[gen_mom]]) == 6)))
        lep_match_dR = deltaR(lep_eta,lep_phi,gen_eta[lep_cut],gen_phi[lep_cut])
        df_[key_]['val']['matchedGenLep'] = (lep_match_dR < .1).sum() > 0 
        #
        if ('TTZH' in key_):
                isbb_fromZ = (((gen_ids == -5) | (gen_ids == 5)) & 
                              (gen_ids[gen_mom] == 23))
                isqq_fromZ = ((abs(gen_ids)<5) & 
                              (gen_ids[gen_mom] == 23))
                isbb_fromH = ((abs(gen_ids) == 5) &
                              (gen_ids[gen_mom] == 25))
                isHbb  = ((gen_ids == 25) & (isbb_fromH.sum() == 2))
                isZbb  = ((gen_ids == 23) & (isbb_fromZ.sum() == 2))
                isZqq  = ((gen_ids == 23) & (isqq_fromZ.sum() == 2))
                isZH = (isHbb | isZbb | isZqq)
                ##### testing things t pt vs z/h pt
                print('isZbb',sum(weight[isZbb.sum() > 0] ))
                print('isZqq',sum(weight[isZqq.sum() > 0] ))
                print('isHbb',sum(weight[isHbb.sum() > 0] ))
                #ist    = (abs(gen_ids) == 6)
                #ist_Zbb= (ist & (isbb_fromZ.sum() > 0))
                #ist_Hbb= (ist & (isHbb.sum() > 0))
                #t_ZbbMaxpt =  np.nanmax(fillne(gen_pt[ist_Zbb]), axis=1)
                #t_ZbbMinpt =  np.nanmin(fillne(gen_pt[ist_Zbb]), axis=1)
                #t_HbbMaxpt =  np.nanmax(fillne(gen_pt[ist_Hbb]), axis=1)
                #t_HbbMinpt =  np.nanmin(fillne(gen_pt[ist_Hbb]), axis=1)
                #
                #Hbbpt = fillne(gen_pt[isHbb]).flatten()
                #Zbbpt = fillne(gen_pt[isZbb]).flatten()
                #print(Hbbpt.shape)
                #print(t_HbbMaxpt.shape)
                #plt.hist2d(x=Hbbpt[Hbbpt >= 200], y=t_HbbMaxpt[Hbbpt >= 200], range = ((200,500),(0,500)), cmin= 0.1, bins= 15)
                #plt.title('t max pt vs H pt')
                #plt.xlabel('H pt')
                #plt.ylabel('t max pt')
                #plt.show()
                #plt.close()
                #plt.hist2d(x=Hbbpt[Hbbpt >= 200], y=t_HbbMinpt[Hbbpt >= 200], range = ((200,500),(0,500)), cmin= 0.1, bins= 15)
                #plt.title('t min pt vs H pt')
                #plt.xlabel('H pt')
                #plt.ylabel('t min pt')
                #plt.show()
                #plt.close()
                #plt.hist2d(x=Zbbpt[Zbbpt >= 200], y=t_ZbbMaxpt[Zbbpt >= 200], range = ((200,500),(0,500)), cmin= 0.1, bins= 15)
                #plt.title('t max pt vs Z pt')
                #plt.xlabel('Z pt')
                #plt.ylabel('t max pt')
                #plt.show()
                #plt.close()
                #plt.hist2d(x=Zbbpt[Zbbpt >= 200], y=t_ZbbMinpt[Zbbpt >= 200], range = ((200,500),(0,500)), cmin= 0.1, bins= 15)
                #plt.title('t min pt vs Z pt')
                #plt.xlabel('Z pt')
                #plt.ylabel('t min pt')
                #plt.show()
                #plt.close()
                #exit()
                #
                pt_cut = int(cfg.skim_ZHbb_dir.split('_')[-1][:3])
                zh_pt = fill1e(gen_pt[isZH]).flatten()
                zh_eta = fill1e(gen_eta[isZH]).flatten()
                zh_phi = fill1e(gen_phi[isZH]).flatten()
                #
                zh_match_dR = deltaR(zh_eta,zh_phi,fj_eta,fj_phi)
                print('isZbb, ZH_pt >= 100',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.))] ))
                print('isZqq, ZH_pt >= 100',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.))] ))
                print('isHbb, ZH_pt >= 100',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.))] ))
                print('isZbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6)] ))
                print('isZqq, ZH_pt >= 100, abs(ZH_eta) <= 2.6',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6)] ))
                print('isHbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6)] ))
                print('isZbb, ZH_pt >= 100, abs(ZH_eta) 2.0 - 2.6',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (abs(zh_eta) >= 2.0)] ))
                print('isZqq, ZH_pt >= 100, abs(ZH_eta) 2.0 - 2.6',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (abs(zh_eta) >= 2.0)] ))
                print('isHbb, ZH_pt >= 100, abs(ZH_eta) 2.0 - 2.6',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (abs(zh_eta) >= 2.0)] ))
                print('isZbb, ZH_pt >= 100, abs(ZH_eta) 2.0 - 2.6, reco matched',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (abs(zh_eta) >= 2.0) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                print('isZqq, ZH_pt >= 100, abs(ZH_eta) 2.0 - 2.6, reco matched',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (abs(zh_eta) >= 2.0) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                print('isHbb, ZH_pt >= 100, abs(ZH_eta) 2.0 - 2.6, reco matched',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (abs(zh_eta) >= 2.0) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                print('isZbb, ZH_pt >= 100, abs(ZH_eta) 0.0 - 2.0',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.0) & (abs(zh_eta) >= 0.0)] ))
                print('isZqq, ZH_pt >= 100, abs(ZH_eta) 0.0 - 2.0',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.0) & (abs(zh_eta) >= 0.0)] ))
                print('isHbb, ZH_pt >= 100, abs(ZH_eta) 0.0 - 2.0',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.0) & (abs(zh_eta) >= 0.0)] ))
                print('isZbb, ZH_pt >= 100, abs(ZH_eta) 0.0 - 2.0, reco matched',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.0) & (abs(zh_eta) >= 0.0) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                print('isZqq, ZH_pt >= 100, abs(ZH_eta) 0.0 - 2.0, reco matched',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.0) & (abs(zh_eta) >= 0.0) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                print('isHbb, ZH_pt >= 100, abs(ZH_eta) 0.0 - 2.0, reco matched',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.0) & (abs(zh_eta) >= 0.0) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                #print('isZbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj,',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & ((zh_match_dR <= 0.8).sum() > 0)] ))
                #print('isZqq, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj,',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & ((zh_match_dR <= 0.8).sum() > 0)] ))
                #print('isHbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj,',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & ((zh_match_dR <= 0.8).sum() > 0)] ))
                #print('isZbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200,',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut)).sum() > 0)] ))
                #print('isZqq, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200,',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut)).sum() > 0)] ))
                #print('isHbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200,',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut)).sum() > 0)] ))
                #print('isZbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200, Hscore >= 0.0,',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                #print('isZqq, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200, Hscore >= 0.0,',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                #print('isHbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200, Hscore >= 0.0,',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0)] ))
                #print('isZbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200, Hscore >= 0.0,lep_match_dr,',sum(weight[(isZbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0) & ((lep_match_dR <= 0.1).sum() > 0)] ))
                #print('isZqq, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200, Hscore >= 0.0,lep_match_dr,',sum(weight[(isZqq.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0) & ((lep_match_dR <= 0.1).sum() > 0)] ))
                #print('isHbb, ZH_pt >= 100, abs(ZH_eta) <= 2.6, dr_matched_fj, fj_pt >= 200, Hscore >= 0.0,lep_match_dr,',sum(weight[(isHbb.sum() > 0) & (zh_pt >= (pt_cut - 100.)) & (abs(zh_eta) <= 2.6) & (((zh_match_dR <= 0.8) & (fj_pt >= pt_cut) & (fj_Hscore >= 0.0)).sum() > 0) & ((lep_match_dR <= 0.1).sum() > 0)] ))
                ## Test
                ind_zh = np.argsort(fillne(zh_match_dR),axis=1)
                zh_match_dr_sort = np.take_along_axis(fillne(zh_match_dR),ind_zh,axis=1)
                fj_pt_sort =  np.take_along_axis(fillne(fj_pt),ind_zh,axis=1)
                dPt = (zh_pt - fj_pt_sort.T).T
                #
                fjpt300_only = np.where((fj_pt_sort>=300),zh_match_dr_sort,np.nan)
                fjpt300_sec = np.where(np.isnan(fjpt300_only[:,0]), fjpt300_only[:,1],np.nan)
                #print(fjpt300_sec)
                #plt.hist(zh_match_dr_sort[:,0], bins = 50, range = (0,5))
                #plt.title('dr of closest fj to gen Z/H')
                #plt.show()
                #plt.close()
                #plt.hist(dPt[:,0], bins=100, range =(-500,500))
                #plt.title('dPt of closest fj to gen Z/H')
                #plt.show()
                #plt.close()
                ##
                #plt.hist(np.where(zh_match_dr_sort < 0.8,dPt, np.nan)[:,0], bins = 100, range = (-500,500))
                #plt.title('dPt of closest fj to gen Z/H (dr < 0.8)')
                #plt.show()
                #plt.close()
                #plt.hist(np.where((zh_match_dr_sort < 0.8) & (fj_pt_sort >= 300),dPt, np.nan)[:,0], bins = 100, range = (-500,500))
                #plt.title('dPt of closest fj (pt >= 300) to gen Z/H (dr < 0.8)')
                #plt.show()
                #plt.close()
                #plt.hist2d(x=np.where((zh_match_dr_sort < 0.8) & (fj_pt_sort >= 300),fj_pt_sort, np.nan)[:,0],y=zh_pt,bins=50,range = ((250,500),(0,500)), cmin=0.1 )
                #plt.show()
                #plt.close()
                ##
                #plt.hist(fj_pt_sort[:,0], bins=50, range= (200,500))
                #plt.title('fj pt of closest fj to gen Z/H')
                #plt.show()
                #plt.close()
                #plt.hist( fjpt300_only[:,0], bins = 50, range = (0,5))
                #plt.title('dr of closest fj (pt >= 300) to gen Z/H')
                #plt.show()
                #plt.close()
                #plt.hist(fjpt300_sec, bins = 50, range = (0,5))
                #plt.title('dr of second closest fj (pt >= 300), where closest fj (pt < 300)')
                #plt.show()
                #plt.close()
                ##
                zh_match = ((fj_pt >= pt_cut ) & (fj_Hscore >= 0.0) & 
                            (zh_match_dR <= 0.8) & (zh_pt >= (pt_cut-100.)) & (zh_eta <= 2.6) & (zh_eta >= -2.6))
                #
                df_[key_]['val']['Zbb']= (isZbb.sum() > 0)
                df_[key_]['val']['Hbb']= (isHbb.sum() > 0)
                df_[key_]['val']['Zqq']= (isZqq.sum() > 0)
                df_[key_]['val']['genZHpt']  = zh_pt
                df_[key_]['val']['genZHeta'] = zh_eta
                #
                df_[key_]['val']['matchedGenZH'] = (zh_match).sum() > 0 
                df_[key_]['val']['matchedGen_ZHbb'] = (((zh_match).sum() > 0) & ((lep_match_dR <= .1).sum() > 0)& (isZqq.sum() == 0))
                df_[key_]['val']['matchedGen_Zbb'] = (((zh_match).sum() > 0) & ((lep_match_dR <= .1).sum() > 0) & (isZbb.sum() > 0))
                df_[key_]['val']['matchedGen_Hbb'] = (((zh_match).sum() > 0) & ((lep_match_dR <= .1).sum() > 0) & (isHbb.sum() > 0))
                df_[key_]['val']['matchedGen_Zqq'] = (((zh_match).sum() > 0) & ((lep_match_dR <= .1).sum() > 0) & (isZqq.sum() > 0))
                print('Matched GenZH',         len(df_[key_]['val']['matchedGenZH']))
                print('Matched GenZHtobb',     sum(df_[key_]['val']['matchedGen_ZHbb']))
                print('Matched GenZtobb',      sum(df_[key_]['val']['matchedGen_Zbb']))
                print('Matched GenHtobb',      sum(df_[key_]['val']['matchedGen_Hbb']))
                print('Matched GenZtoqq',      sum(df_[key_]['val']['matchedGen_Zqq']))
                print('No matched GenZHtobb',  sum((df_[key_]['val']['matchedGen_ZHbb'] == False) & (df_[key_]['val']['matchedGen_Zqq'] == False)))
                #
                print('# w_events GenZH',     sum(weight[df_[key_]['val']['matchedGenZH']]))
                print('# w_events GenZHtobb', sum(weight[df_[key_]['val']['matchedGen_ZHbb']]))
                print('# w_events GenZtobb',  sum(weight[df_[key_]['val']['matchedGen_Zbb']]))
                print('# w_events GenHtobb',  sum(weight[df_[key_]['val']['matchedGen_Hbb']]))
                print('# w_events GenZtoqq',  sum(weight[df_[key_]['val']['matchedGen_Zqq']]))
                print('# w_events GenZHtobb', sum(weight[(df_[key_]['val']['matchedGen_ZHbb'] == False) & (df_[key_]['val']['matchedGen_Zqq'] == False) & (isZqq.sum() == 0)]))
        #
        df_[key_]['val'].to_pickle(outDir_+'result_'+year+'_'+sample+'_val.pkl')

def lepCleaned_v2(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df_ = kFit.retrieveData(files_, samples_, outDir_, getak8_=True)
    for key_ in df_.keys():
        df = df_[key_]
        fat_df = df['ak8']
        jet_df = df['df']
        val_df = df['val']
        #
        lep_eta  = val_df['Lep_eta'].values
        lep_phi  = val_df['Lep_phi'].values
        #
        fj_eta = fat_df['eta']
        fj_phi = fat_df['phi']
        #
        _eta = []
        _phi = []
        for key_str in jet_df.keys():
            if   ('eta_' in key_str):
                _eta.append(key_str)
            elif ('phi_' in key_str):
                _phi.append(key_str)
        j_eta = jet_df[_eta].to_numpy()
        j_phi = jet_df[_phi].to_numpy()
        ##
        # Clean Fat Jets
        ##
        lep_fj_dR = deltaR(lep_eta,lep_phi,fj_eta, fj_phi)
        df_[key_]['ak8']['fj_lep_mask'] = lep_fj_dR > .8
        ##
        # Clean ak4 Jets
        ## 
        lep_j_dR = deltaR(lep_eta,lep_phi,j_eta, j_phi)
        df_[key_]['val']['j_lep_mask'] = (lep_j_dR > .4).tolist()
        # Save file
        year = key_.split('_')[1]
        sample = key_.split('_')[0]
        with open(outDir_+'result_'+year+'_'+sample+'_ak8.pkl', 'wb') as handle:
            pickle.dump(df_[key_]['ak8'], handle, protocol=pickle.HIGHEST_PROTOCOL)
        df_[key_]['val'].to_pickle(outDir_+'result_'+year+'_'+sample+'_val.pkl')
        
        
def fixttH_weight(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap):
    df_ = kFit.retrieveData(files_, ['TTZH'], outDir_, getgen_=True)
    df = df_['TTZH_2017']
    w  = df['val']['weight'].values
    
    #
    gen_df = df['gen']
    zh_ids = gen_df['GenPart_pdgId']
    isHbb = ((zh_ids == 25).sum() > 0)
    #
    # divide by current bad cross-section (.153) and multiply by the correct cross section (.2934)
    # original weight stored, this is to keep from double weighting if run over more than 2x
    print(w[isHbb])
    fixed_w = np.where(isHbb,np.sign(w)*0.001127809*(.2934/.153),w)
    df_['TTZH_2017']['val']['weight'] = fixed_w
    df_['TTZH_2017']['val'].to_pickle(outDir_+'result_2017_TTZH_val.pkl')
    
if __name__ == '__main__':   

    files_samples_outDir = cfg.ZHbbFitCfg
    #
    #prD.getData(         *files_samples_outDir, *cfg.ZHbbFitCut, cfg.ZHbbFitMaxJets, treeDir_ = cfg.tree_dir+'_bb', getGenData_ = True, getak8var_=True)
    ####prD.interpData(      *files_samples_outDir, cfg.ZHbbFitMaxJets)  
    #fixttH_weight(*files_samples_outDir, cfg.ZHbbFitoverlap)
    #lepCleaned_v2(*files_samples_outDir, cfg.ZHbbFitoverlap)
    #matchLep(*files_samples_outDir, cfg.ZHbbFitoverlap)
    ###
    #ZHbbAna(*files_samples_outDir, cfg.ZHbbFitoverlap)
    ########
    #GenAna(*files_samples_outDir, cfg.ZHbbFitoverlap)
    #GenAna_ttbar(*files_samples_outDir, cfg.ZHbbFitoverlap)
    ########
    plotAna(*files_samples_outDir, cfg.ZHbbFitoverlap)
    #Bkg_Est(*files_samples_outDir, cfg.ZHbbFitoverlap)
