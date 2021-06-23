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

@t2Run
def reco_soft_lep_helper(obj):
    '''
    use this to calculate the min(invariant M (sel lep, sfos soft lep))
    '''
    getTLVm  = TLorentzVectorArray.from_ptetaphim  
    common_k = ['pt','eta','phi','mass']
    # Electrons 
    sel_elec = [obj.val_df[obj.val_df['passSingleLepElec'] == 1]['Lep_'+k] for k in common_k]
    sel_elec_ch = obj.val_df[obj.val_df['passSingleLepElec'] == 1]['Lep_ch']
    sel_elec_tlv = getTLVm(*sel_elec)
    soft_elec = [obj.softe_df[obj.val_df['passSingleLepElec'] == 1]['Electron_'+k] for k in common_k]
    soft_elec_ch = obj.softe_df[obj.val_df['passSingleLepElec'] == 1]['Electron_charge']
    soft_elec_tlv = getTLVm(*soft_elec)

    #print(sum((~(soft_elec_tlv == sel_elec_tlv)).counts) > 0)
    if sum((~(soft_elec_tlv == sel_elec_tlv)).counts) > 0 :
        sel_soft_elec_invm = (sel_elec_tlv+soft_elec_tlv[(~(soft_elec_tlv == sel_elec_tlv) & ~(sel_elec_ch == soft_elec_ch))]).mass
    else:
        sel_soft_elec_invm = np.array([9939])
    #sel_soft_elec_invm = (sel_elec_tlv+soft_elec_tlv[~(soft_elec_tlv == sel_elec_tlv)]).mass
        # Muons
    sel_mu = [obj.val_df[obj.val_df['passSingleLepMu'] == 1]['Lep_'+k] for k in common_k]
    sel_mu_ch = obj.val_df[obj.val_df['passSingleLepMu'] == 1]['Lep_ch']
    sel_mu_tlv = getTLVm(*sel_mu)
    soft_mu = [obj.softmu_df[obj.val_df['passSingleLepMu'] == 1]['Muon_'+k] for k in common_k]
    soft_mu_ch = obj.softmu_df[obj.val_df['passSingleLepMu'] == 1]['Muon_charge']
    soft_mu_tlv = getTLVm(*soft_mu)
    #
    if sum((~(soft_elec_tlv == sel_elec_tlv)).counts) > 0 :
        sel_soft_mu_invm = (sel_mu_tlv+soft_mu_tlv[(~(soft_mu_tlv == sel_mu_tlv) & ~(sel_mu_ch == soft_mu_ch))]).mass
    else:
        sel_soft_mu_invm = np.array([999])
    #sel_soft_mu_invm = (sel_mu_tlv+soft_mu_tlv[~(soft_mu_tlv == sel_mu_tlv)]).mass
    # 
    #obj.val_df['min_sel_soft_elec_invm'] = np.where(obj.val_df['passSingleLepElec'] == 1, 0, 999)
    obj.val_df['min_sel_soft_elec_invm'] = 999
    obj.val_df.loc[obj.val_df['passSingleLepElec'] == 1,['min_sel_soft_elec_invm']] = np.where(np.isinf(sel_soft_elec_invm.min()), 999, sel_soft_elec_invm.min())
    #obj.val_df['min_sel_soft_mu_invm'] = np.where(obj.val_df['passSingleLepMu'] == 1, 0, 999)
    obj.val_df['min_sel_soft_mu_invm'] = 999
    obj.val_df.loc[obj.val_df['passSingleLepMu'] == 1,['min_sel_soft_mu_invm']]   = np.where(np.isinf(sel_soft_mu_invm.min()), 999, sel_soft_mu_invm.min())
    # Combine to one
    obj.val_df['passNotHadLep'] = False
    obj.val_df.loc[obj.val_df['passSingleLepElec'] == 1, ['passNotHadLep']] =  (obj.val_df['min_sel_soft_elec_invm'][obj.val_df['passSingleLepElec'] == 1] > 12)
    obj.val_df.loc[obj.val_df['passSingleLepMu']   == 1, ['passNotHadLep']] =  (obj.val_df['min_sel_soft_mu_invm'][  obj.val_df['passSingleLepMu']   == 1] > 12)
    
