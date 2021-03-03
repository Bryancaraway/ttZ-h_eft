import os
import sys
from uproot_methods import TLorentzVectorArray 
import config.ana_cff as cfg
import lib.fun_library as lib
from lib.fun_library import fill1e, fillne, deltaR, deltaPhi, invM, calc_mtb, t2Run
from modules.AnaDict import AnaDict
#
import numpy as np
import pandas as pd


def reco_zh_helper_andrew(obj):
    # Note: in order to take advatage of numpy methods with dimensionallity
    # it is necessary to numpy-ize (rectangularize) jagged-arrays using custom (hacky)
    # function 'fillne'
    k_getter = (lambda df, keys: df.loc(keys).values() )
    # create event level variables here
    # get ak4, ak8,val variables needed for analysis #
    l = 'Lep_'
    j = 'Jet_'
    fj= 'FatJet_'
    common_k = ['pt','eta','phi']
    #ak4_pt, ak4_eta, ak4_phi, ak4_mass, ak4_btag = obj.ak4_df.loc(
    #    [j+str_+obj.lc for str_ in ['pt','eta','phi','mass','btagDeepB']]
    #)[obj.ak4_df[j+'lep_mask']].values()
    ak4_pt, ak4_eta, ak4_phi, ak4_mass, ak4_btag = k_getter(
     obj.ak4_df, [j+str_ for str_ in common_k+['mass','btagDeepB']]
    )
    ht= ak4_pt.sum()
    #print(ht)
    obj.val_df['HT'] = ht
    b_pt, b_eta, b_phi, b_mass, b_btag = [ak4_k[ak4_btag >= obj.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_mass,ak4_btag]]
    q_pt, q_eta, q_phi, q_mass, q_btag = [ak4_k[ak4_btag <  obj.b_wp] for ak4_k in [ak4_pt,ak4_eta,ak4_phi,ak4_mass,ak4_btag]]
    #ak8_pt, ak8_eta, ak8_phi, sd_M, ak8_bbtag, ak8_Zhbbtag, w_tag, t_tag, subj1, subj2 = obj.ak8_df.loc(
    #    [fj+str_+obj.lc for str_ in ['pt','eta','phi','msoftdrop','btagDeepB','deepTagMD_bbvsLight','deepTagMD_WvsQCD','deepTagMD_TvsQCD','subJetIdx1','subJetIdx2']]
    #)[obj.ak8_df[fj+'lep_mask']].values()
    ak8_pt, ak8_eta, ak8_phi, sd_M, ak8_bbtag, ak8_Zhbbtag, ak8_doubleB, w_tag, t_tag, subj1, subj2 = k_getter( 
        obj.ak8_df, [fj+str_ for str_ in common_k+['msoftdrop','btagDeepB','deepTagMD_bbvsLight', 'btagHbb', 'deepTagMD_WvsQCD','deepTagMD_TvsQCD','subJetIdx1','subJetIdx2']]
    )
    #subj_pt, subj_eta, subj_phi, subj_btag = obj.ak8_df['SubJet_pt'], obj.ak8_df['SubJet_eta'], obj.ak8_df['SubJet_phi'], obj.ak8_df['SubJet_btagDeepB']
    subj_pt,  subj_btag = obj.ak8_df['SubJet_pt'],  obj.ak8_df['SubJet_btagDeepB']
    #
    lep_pt, lep_eta, lep_phi, lep_M = [obj.val_df[l+str_] for str_ in common_k+['mass']]
    met_pt, met_phi = obj.val_df['MET_pt'], obj.val_df['MET_phi']
    # reconstruct z, h -> bb candidate
    Zh_reco_cut = ((ak8_pt >= obj.pt_cut) & (sd_M >= 50) & (sd_M <= 200) & (ak8_Zhbbtag >= 0.0)) # in future, put kinem cut values in cfg file
    #print(Zh_reco_cut)
    Zh_Zhbbtag,Zh_pt,Zh_eta,Zh_phi,Zh_M, Zh_doubleB, Zh_wtag, Zh_ttag,Zh_bbtag=lib.sortbyscore(
        [ak8_Zhbbtag,ak8_pt,ak8_eta,ak8_phi,sd_M, ak8_doubleB, w_tag,t_tag,ak8_bbtag],ak8_Zhbbtag,Zh_reco_cut)
    # =============================== # 
    # compute subject info related to Zh candidate : have to sort out candidate with no subjet
    #ak8_sj1_pt, ak8_sj2_pt, ak8_sj1_eta, ak8_sj2_eta, ak8_sj1_phi, ak8_sj2_phi, ak8_sj1_btag, ak8_sj2_btag = [subj_k[id_[id_ != -1]] for subj_k in [subj_pt, subj_eta, subj_phi, subj_btag] for id_ in [subj1, subj2]] # sj kinems per ak8jet
    ak8_sj1_pt, ak8_sj2_pt, ak8_sj1_btag, ak8_sj2_btag = [subj_k[id_[id_ != -1]] for subj_k in [subj_pt,  subj_btag] for id_ in [subj1, subj2]] # sj kinems per ak8jet
    ak8_sj1_sdM, ak8_sj2_sdM, ak8_sj1_Zhbb, ak8_sj2_Zhbb = [ak8_k[id_ != -1] for ak8_k in [sd_M, ak8_Zhbbtag] for id_ in [subj1, subj2]] # ak8 kinems 
    #Zh_kinem_sj1_cut = ((ak8_sj1_pt>=obj.pt_cut) & (ak8_sj1_sdM >= 50) & (ak8_sj1_sdM <= 200) & (ak8_sj1_Zhbb >= 0.0))

    Zh_kinem_sj1_cut = (ak8_sj1_Zhbb >= 0.0)
    Zh_sj1_pt, Zh_sj1_btag, Zh_sj1_Zhbb = lib.sortbyscore([
        ak8_sj1_pt,ak8_sj1_btag,ak8_sj1_Zhbb],ak8_sj1_Zhbb,Zh_kinem_sj1_cut)
    #Zh_kinem_sj2_cut = ((ak8_sj2_pt>=obj.pt_cut) & (ak8_sj2_sdM >= 50) & (ak8_sj2_sdM <= 200) & (ak8_sj2_Zhbb >= 0.0))
    Zh_kinem_sj2_cut = ((ak8_sj2_Zhbb >= 0.0))
    Zh_sj2_pt,  Zh_sj2_btag, Zh_sj2_Zhbb = lib.sortbyscore([
        ak8_sj2_pt,ak8_sj2_btag,ak8_sj2_Zhbb],ak8_sj2_Zhbb,Zh_kinem_sj2_cut)
    #
    Zh_sj_b12  = np.column_stack([Zh_sj1_btag[:,0], Zh_sj2_btag[:,0]])
    Zh_sj_pt12 = np.nan_to_num(np.column_stack([Zh_sj1_pt[:,0], Zh_sj2_pt[:,0]]) )

    #Zh_sj1_sj2_dr = deltaR(Zh_sj1_eta, Zh_sj1_phi, Zh_sj2_eta, Zh_sj2_phi)
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
    Zh_b_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(b_eta),fillne(b_phi))
    process_sorted = (lambda k_,dr_,ind_: np.take_along_axis(np.where(dr_ > 0.8, fillne(k_),np.nan),ind_,axis=1) )
    ind_Zh_b_dr = np.argsort(np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis=1)
    b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort  = [
        process_sorted(b_k, Zh_b_dr, ind_Zh_b_dr) for b_k in [b_pt,b_eta,b_phi,b_mass]]
    ind_Zh_b_pt = np.argsort(-1*np.where(Zh_b_dr > 0.8, fillne(b_pt), np.nan),axis=1)
    b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort, b_btag_ptsort  = [
        process_sorted(b_k, Zh_b_dr, ind_Zh_b_pt) for b_k in [b_pt,b_eta,b_phi,b_mass, b_btag]]
    #
    Zh_q_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],fillne(q_eta),fillne(q_phi))
    ind_Zh_q_dr = np.argsort(np.where(Zh_q_dr > 0.8, Zh_q_dr, np.nan),axis=1)
    q_pt_dRsort, q_eta_dRsort, q_phi_dRsort, q_mass_dRsort = [
        process_sorted(q_k, Zh_q_dr, ind_Zh_q_dr) for q_k in [q_pt,q_eta,q_phi,q_mass]]
    ind_Zh_q_pt = np.argsort(-1*np.where(Zh_q_dr > 0.8, fillne(q_pt), np.nan),axis=1)
    q_pt_ptsort, q_eta_ptsort, q_phi_ptsort, q_mass_ptsort, q_btag_ptsort  = [
        process_sorted(q_k, Zh_q_dr, ind_Zh_q_pt) for q_k in [q_pt,q_eta,q_phi,q_mass, q_btag]]

    Zh_l_dr = deltaR(Zh_eta[:,0],Zh_phi[:,0],lep_eta,lep_phi)
    Zh_l_invM_sd = lib.invM(Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],lep_pt,lep_eta,lep_phi,lep_M)
    #l_b_dr = deltaR(lep_eta.values,lep_phi.values,*map(fillne,[b_eta,b_phi]))
    b_b_dr = deltaR(b_eta_ptsort[:,0],b_phi_ptsort[:,0], b_eta_ptsort[:,1],b_phi_ptsort[:,1])
    q_q_dr = deltaR(q_eta_ptsort[:,0],q_phi_ptsort[:,0], q_eta_ptsort[:,1],q_phi_ptsort[:,1])
    l_b_dr = deltaR(lep_eta.values,lep_phi.values, b_eta_ptsort,b_phi_ptsort)
    #l_b_invM = lib.invM(lep_pt.values,lep_eta.values,lep_phi.values,lep_M.values,*map(fillne,[b_pt,b_eta,b_phi,b_mass]))
    l_b_invM = lib.invM(lep_pt.values,lep_eta.values,lep_phi.values,lep_M.values, b_pt_ptsort,b_eta_ptsort,b_phi_ptsort,b_mass_ptsort)
    #
    getTLV  = TLorentzVectorArray.from_ptetaphi
    getTLVm = TLorentzVectorArray.from_ptetaphim
    #
    #b_tlv   = getTLVm(*map(lambda p: p.T,[b_pt_ptsort,b_eta_ptsort,b_phi_ptsort,b_mass_ptsort]) )
    b1_tlv = getTLVm(b_pt_ptsort[:,0],b_eta_ptsort[:,0],b_phi_ptsort[:,0],b_mass_ptsort[:,0])
    b2_tlv = getTLVm(b_pt_ptsort[:,1],b_eta_ptsort[:,1],b_phi_ptsort[:,1],b_mass_ptsort[:,1])
    q_tlv   = getTLVm(*map(lambda p: p.T,[q_pt_ptsort,q_eta_ptsort,q_phi_ptsort,q_mass_ptsort]) )
    lep_tlv = getTLV( lep_pt.values,lep_eta.values,lep_phi.values,lep_M.values)
    #
    ind_lb = np.argsort(l_b_dr,axis=1) 
    l_b_invM_dRsort, l_b_dr_dRsort = [np.take_along_axis(lb_comb,ind_lb,axis=1) for lb_comb in [l_b_invM,l_b_dr]]
    max_l_b_dr,  min_l_b_dr  = np.nanmax(l_b_dr_dRsort, axis=1),     np.nanmin(l_b_dr_dRsort, axis=1)
    max_lb_invM, min_lb_invM = np.nanmax(l_b_invM_dRsort, axis=1), np.nanmin(l_b_invM_dRsort, axis=1)
    # farthest b to l
    for i, btlv in enumerate([b1_tlv,b2_tlv]):
        obj.val_df[f'outZH_b{i+1}_pt']    = b_pt_ptsort[:,i]
        obj.val_df[f'outZH_b{i+1}_score'] = b_btag_ptsort[:,i]
        b_l_tlv = lep_tlv + btlv
        obj.val_df[f'l_b{i+1}_mtb'] = lib.calc_mtb(b_l_tlv.pt,b_l_tlv.phi,met_pt.values,met_phi.values)            
        b_q_dr = deltaR(b_eta_ptsort[:,i],b_phi_ptsort[:,i], q_eta_ptsort,q_phi_ptsort) # get 1st and 2nd closest quarks
        obj.val_df[f'outZH_b{i+1}_q_mindr'] = np.nan_to_num(np.nanmin(b_q_dr, axis=1)) # new
        ind_bq_dr = np.argsort(b_q_dr,axis=1)
        q_pt_bdRsort, q_eta_bdRsort, q_phi_bdRsort, q_mass_bdRsort = [np.take_along_axis(q_,ind_bq_dr,axis=1) for q_ in [q_pt_ptsort,q_eta_ptsort,q_phi_ptsort,q_mass_ptsort]]
        q1_tlv = getTLVm(q_pt_bdRsort[:,0],q_eta_bdRsort[:,0],q_phi_bdRsort[:,0],q_mass_bdRsort[:,0])
        q2_tlv = getTLVm(q_pt_bdRsort[:,1],q_eta_bdRsort[:,1],q_phi_bdRsort[:,1],q_mass_bdRsort[:,1])
        obj.val_df[f'outZH_q_q_dr_nearb{i+1}'] = np.nan_to_num(q1_tlv.delta_r(q2_tlv))
        obj.val_df[f'outZH_qq_M_nearb{i+1}'] = np.nan_to_num((q1_tlv + q2_tlv).mass)
        obj.val_df[f'outZH_b{i+1}q_M'] = np.nan_to_num((btlv + q1_tlv).mass)
        bqq_tlv = btlv + q1_tlv + q2_tlv
        obj.val_df[f'outZH_b{i+1}_qq_dr'] = np.nan_to_num(btlv.delta_r(q1_tlv + q2_tlv))
        obj.val_df[f'outZH_b{i+1}qq_M'] = np.nan_to_num(bqq_tlv.mass) # new
        lbbqq_tlv = lep_tlv + btlv + b2_tlv + q1_tlv + q2_tlv
        obj.val_df[f'ZH_b{i+1}qq_dr'] = np.nan_to_num(deltaR(Zh_eta[:,0],Zh_phi[:,0],bqq_tlv.eta, bqq_tlv.phi)) # new
        obj.val_df[f'ZH_lbb{i+1}qq_dr'] = np.nan_to_num(deltaR(Zh_eta[:,0],Zh_phi[:,0],lbbqq_tlv.eta, lbbqq_tlv.phi)) # new
        #
    #
    n_b_outZhbb = np.nansum(Zh_b_dr > .8, axis=1)
    n_q_outZhbb = np.nansum(Zh_q_dr > .8 , axis=1)
    n_b_inZhbb  = np.nansum(Zh_b_dr <= .8, axis=1)
    n_q_inZhbb  = np.nansum(Zh_q_dr <= .8, axis=1)
    
    ind_inZh_b_dr = np.argsort(np.where(Zh_b_dr <= 0.8, Zh_b_dr, np.nan),axis=1)
    b_pt_indRsort, b_eta_indRsort, b_phi_indRsort, b_mass_indRsort  = [
        np.take_along_axis(np.where(Zh_b_dr <= 0.8, fillne(b_k),np.nan),ind_inZh_b_dr,axis=1) for b_k in [b_pt,b_eta,b_phi,b_mass]]
    inZhb_outZhb_dr = np.nan_to_num(deltaR(b_eta_indRsort[:,0], b_phi_indRsort[:,0], b_eta_dRsort[:,0], b_phi_dRsort[:,0]))

    #Zh_ak4_dr, Zh_b_dr, Zh_q_dr = fillne(Zh_ak4_dr), fillne(Zh_b_dr), fillne(Zh_q_dr)
    b12_pt_ptsort =  (b1_tlv+b2_tlv).pt
    b12_m_ptsort  =  (b1_tlv+b2_tlv).mass # new
    b12_dr_ptsort =  b1_tlv.delta_r(b2_tlv) # new 
    h_t_b         = np.nansum(b_pt_ptsort, axis=1) # new
    #
    sc_pt_outZh   = h_t_b + np.nansum(q_pt_ptsort, axis=1) + lep_pt # new
    #
    obj.val_df['Zh_closeb_invM'] = lib.invM(
        Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],b_pt_dRsort[:,0],b_eta_dRsort[:,0],b_phi_dRsort[:,0],b_mass_dRsort[:,0])
    obj.val_df['Zh_2ndcloseb_invM'] = lib.invM(
        Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],b_pt_dRsort[:,1],b_eta_dRsort[:,1],b_phi_dRsort[:,1],b_mass_dRsort[:,1])
    obj.val_df['Zh_closeq_invM'] = np.nan_to_num(
        lib.invM(
            Zh_pt[:,0],Zh_eta[:,0],Zh_phi[:,0], Zh_M[:,0],q_pt_dRsort[:,0],q_eta_dRsort[:,0],q_phi_dRsort[:,0],q_mass_dRsort[:,0]))
    #
    ind_Zh_b_dr_nocut = np.argsort(Zh_b_dr,axis=1)
    b_pt_ncdRsort, b_eta_ncdRsort, b_phi_ncdRsort, b_mass_ncdRsort  = [
        np.take_along_axis(fillne(b_k), ind_Zh_b_dr_nocut,axis=1) for b_k in [b_pt,b_eta,b_phi,b_mass]]
    ind_Zh_q_dr_nocut = np.argsort(Zh_q_dr,axis=1)
    q_pt_ncdRsort, q_eta_ncdRsort, q_phi_ncdRsort, q_mass_ncdRsort  = [
        np.take_along_axis(fillne(q_k), ind_Zh_q_dr_nocut,axis=1) for q_k in [q_pt,q_eta,q_phi,q_mass]]
    #
    obj.val_df['Zh_closeb_dr'] = deltaR(
        Zh_eta[:,0],Zh_phi[:,0],b_eta_ncdRsort[:,0],b_phi_ncdRsort[:,0])
    obj.val_df['Zh_2ndcloseb_dr'] = deltaR(
        Zh_eta[:,0],Zh_phi[:,0],b_eta_ncdRsort[:,1],b_phi_ncdRsort[:,1])
    obj.val_df['Zh_closeq_dr'] = np.nan_to_num(deltaR(
        Zh_eta[:,0],Zh_phi[:,0],q_eta_ncdRsort[:,0],q_phi_ncdRsort[:,0]))

    #
    obj.val_df['nBottoms']  = b_pt.counts
    obj.val_df['n_ak4jets'] = ak4_pt.counts
    obj.val_df['n_ak8jets'] = ak8_pt.counts
    #
    obj.val_df['jetpt_1']  = ak4_pt.pad(2)[:,0]
    obj.val_df['jetpt_2']  = ak4_pt.pad(2)[:,1]
    obj.val_df['bjetpt_1']  = b_pt.pad(2)[:,0]
    obj.val_df['bjetpt_2']  = b_pt.pad(2)[:,1]
    obj.val_df['jeteta_1'] = ak4_eta.pad(2)[:,0]
    obj.val_df['jeteta_2'] = ak4_eta.pad(2)[:,1]
    obj.val_df['bjeteta_1'] = b_eta.pad(2)[:,0]
    obj.val_df['bjeteta_2'] = b_eta.pad(2)[:,1]
    obj.val_df['jetbtag_1'] = ak4_btag.pad(2)[:,0]
    obj.val_df['jetbtag_2'] = ak4_btag.pad(2)[:,1]
    obj.val_df['bjetbtag_1'] = b_btag.pad(2)[:,0]
    obj.val_df['bjetbtag_2'] = b_btag.pad(2)[:,1]
    #
    obj.val_df['fjetpt_1']  = ak8_pt.pad(1)[:,0]
    obj.val_df['fjeteta_1'] = ak8_eta.pad(1)[:,0]
    obj.val_df['fjetsdm_1'] = sd_M.pad(1)[:,0]
    obj.val_df['fjetbbvl_1']= ak8_Zhbbtag.pad(1)[:,0]
    obj.val_df['fjetwscore_1']= w_tag.pad(1)[:,0]
    obj.val_df['fjettscore_1']= t_tag.pad(1)[:,0]
    #
    obj.val_df['n_ak8_Zhbb'] = ak8_pt[((ak8_pt >= obj.pt_cut) & (sd_M >= 50) & (sd_M <= 200) & (ak8_Zhbbtag >= 0.6))].counts
    obj.val_df['Zh_bbvLscore']  = Zh_Zhbbtag[:,0]
    obj.val_df['Zh_pt']     = Zh_pt[:,0]
    obj.val_df['Zh_eta']    = Zh_eta[:,0]
    obj.val_df['Zh_phi']    = Zh_phi[:,0]
    obj.val_df['Zh_M']      = Zh_M[:,0]
    obj.val_df['Zh_doubleB']      = Zh_doubleB[:,0]
    obj.val_df['Zh_Wscore'] = Zh_wtag[:,0]
    obj.val_df['Zh_Tscore'] = Zh_ttag[:,0]
    obj.val_df['Zh_deepB']  = Zh_bbtag[:,0]
    obj.val_df['n_Zh_btag_sj'] = np.nansum(Zh_sj_b12 >= obj.b_wp, axis=1)
    obj.val_df['n_Zh_sj']       = np.nansum(Zh_sj_b12 >= 0, axis=1)
    obj.val_df['Zh_bestb_sj']   = np.nan_to_num(np.nanmax(Zh_sj_b12, axis=1))
    obj.val_df['Zh_worstb_sj']  = np.nan_to_num(np.nanmin(Zh_sj_b12, axis=1))
    #obj.val_df['Zh_sj1_sj2_dr']    = Zh_sj1_sj2_dr
    obj.val_df['sjpt1_over_Zhpt'] = Zh_sjpt1_over_fjpt
    obj.val_df['sjpt2_over_Zhpt'] = Zh_sjpt2_over_fjpt
    obj.val_df['outZh_max_Wscore'] = np.max(np.nan_to_num(Zh_wtag[:,1:]), axis=1) # new
    obj.val_df['outZh_max_Tscore'] = np.max(np.nan_to_num(Zh_ttag[:,1:]), axis=1) # new
    obj.val_df['outZh_max_bbvLscore']  = np.max(np.nan_to_num(Zh_Zhbbtag[:,1:]), axis=1) # new
    obj.val_df['outZh_max_ak8pt']  = np.max(np.nan_to_num(Zh_pt[:,1:]), axis=1) # new
    obj.val_df['outZh_max_ak8sdM']  = np.max(np.nan_to_num(Zh_M[:,1:]), axis=1) # new

    obj.val_df['spher'] = spher
    obj.val_df['aplan'] = aplan
    
    obj.val_df['max_lb_dr'] = max_l_b_dr
    obj.val_df['min_lb_dr'] = min_l_b_dr
    obj.val_df['max_lb_invM'] = max_lb_invM
    obj.val_df['min_lb_invM'] = min_lb_invM

    obj.val_df['inZhb_outZhb_dr'] = inZhb_outZhb_dr
    obj.val_df['outZh_b12_m'] = b12_m_ptsort
    obj.val_df['outZh_b12_dr'] = b12_dr_ptsort
    obj.val_df['ht_b']         = h_t_b # new
    obj.val_df['ht_outZh']     = sc_pt_outZh # new
    obj.val_df['nonZhbb_q1_dr'] =  np.nan_to_num(np.sort( np.where(Zh_q_dr > 0.8, Zh_q_dr, np.nan),axis = 1 )[:,0])
    obj.val_df['nonZhbb_b1_dr'] =  np.sort( np.where(Zh_b_dr > 0.8, Zh_b_dr, np.nan),axis = 1 )[:,0]
    obj.val_df['n_b_outZh']     = n_b_outZhbb
    obj.val_df['n_b_inZh']      = n_b_inZhbb
    obj.val_df['n_q_outZh']     = n_q_outZhbb
    obj.val_df['n_q_inZh']      = n_q_inZhbb
    obj.val_df['Zh_l_dr']       = Zh_l_dr
    obj.val_df['Zh_l_invM_sd']  = Zh_l_invM_sd
    for i in range(2):
        obj.val_df[f'l_b{i+1}_invM'] = l_b_invM_dRsort[:,i]
        obj.val_df[f'l_b{i+1}_dr']   = l_b_dr_dRsort[:,i]
        #
        obj.val_df[f'outZh_q{i+1}_pt']    = np.nan_to_num(q_pt_ptsort[:,i])
        obj.val_df[f'outZh_q{i+1}_btag']  = np.nan_to_num(q_btag_ptsort[:,i])
