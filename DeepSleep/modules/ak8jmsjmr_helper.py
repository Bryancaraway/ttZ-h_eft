import os
import sys
from uproot_methods import TLorentzVectorArray 
import config.ana_cff as cfg
import lib.fun_library as lib
from lib.fun_library import argmatch, t2Run
from awkward import JaggedArray as aj
from modules.AnaDict import AnaDict
#
import numpy as np
import pandas as pd
import re


# corrections 
#JMS, JMR vals
jmsVals  = {
    '2016':[1.00, 0.9906, 1.0094],
    '2017':[0.982, 0.978, 0.986],
    '2018':[0.982, 0.978, 0.986],
}
jmrVals  = {
    '2016':[1.0, 1.2, 0.8],
    #'2017':[1.09, 1.14, 1.04],
    #'2018':[1.09, 1.14, 1.04],
    '2017':[1.00, 1.05, 0.95],
    '2018':[1.00, 1.05, 0.95],
}
# w-jms
puppicorr_gen = { 
    '2016':(lambda x : 1.0062610-1.0616051*np.power(x*0.079990008,-1.2045377)),
    '2017':(lambda x : 1+0.0090283*np.power(x,(-2*(0.0099852)))-7.30123*np.power(x,-1)),
    '2018':(lambda x : 1+0.0231069*np.power(x,(-2*(0.0704761)))-8.42917*np.power(x,-1))  
}
puppicorr_reco0eta1p3 = {
    '2016':(lambda x: 1.0930198-0.00015006789*x+(3.4486612e-07)*np.power(x,2)+(-2.6810031e-10)*np.power(x,3)+(8.6744023e-14)*np.power(x,4)+(-1.0011359e-17)*np.power(x,5)),
    '2017':(lambda x: 1.04323417805+(8.20581677106e-05)*x+(-2.23790959145e-08)*np.power(x,2)+(-5.56816212196e-12)*np.power(x,3)+(-2.42702058503e-17)*np.power(x,4)+(5.23731618031e-19)*np.power(x,5)),
    '2018':(lambda x: 1.06263222851+(2.97332221436e-05)*x+(-7.31785851974e-09)*np.power(x,2)+(2.53798731754e-13)*np.power(x,3)+(1.68571767997e-16)*np.power(x,4)+(-6.77123709437e-20)*np.power(x,5))
}
puppicorr_reco1p3eta2p5 = {
    '2016':(lambda x: 1.2721152-0.00057164036*x+(8.3728941e-07)*np.power(x,2)+(-5.2043320e-10)*np.power(x,3)+(1.4537521e-13)*np.power(x,4)+(-1.5038869e-17)*np.power(x,5)),
    '2017':(lambda x: 1.11549406241+(-2.01425972518e-05)*x+(8.36181961894e-09)*np.power(x,2)+(4.39451437171e-11)*np.power(x,3)+(1.04302756829e-14)*np.power(x,4)+(-2.10404344784e-17)*np.power(x,5)),
    '2018':(lambda x: 1.11889161475+(2.68579882197e-05)*x+(-4.30234840782e-09)*np.power(x,2)+(8.27000377942e-12)*np.power(x,3)+(1.45823446695e-15)*np.power(x,4)+(-1.65979484436e-18)*np.power(x,5))
}
# w-jmr
puppicorr_massReso_0eta1p3   = (lambda x: 1.0927351+(4.1426227e-05)*x+(-1.3736806e-07)*np.power(x,2)+(1.2295818e-10)*np.power(x,3)+(-4.1970754e-14)*np.power(x,4)+(4.9237927e-18)*np.power(x,5))
puppicorr_massReso_1p3eta2p5 = (lambda x: 1.1649278+(-0.00012678903)*x+(1.0594037e-07)*np.power(x,2)+(6.0720870e-12)*np.power(x,3)+(-1.9924275e-14)*np.power(x,4)+(3.6440065e-18)*np.power(x,5))

@t2Run
def ak8jmsjmr_helper(obj, jec_sys):
    # apply some preliminary cleaning to fatjet variables
    init_cut = (
        (obj.tmp_fatjets["FatJet_pt"]>20) &
        (obj.tmp_fatjets['FatJet_msoftdrop'] >= 20.) & 
        (abs(obj.tmp_fatjets['FatJet_eta']) < 2.4) & 
        ( obj.tmp_fatjets['FatJet_jetId'] >= 2)
    )
    obj.fatjets     = obj.fatjets[init_cut]    
    obj.tmp_fatjets = obj.tmp_fatjets[init_cut]
    # data/mc puppi correction
    jms_helper(obj)
    obj.tmp_fatjets["FatJet_msoftdrop_altnosmear"] = obj.tmp_fatjets["FatJet_msoftdrop"]
    obj.tmp_fatjets["FatJet_msoftdrop_nosmear"] = obj.tmp_fatjets["FatJet_msoftdrop_corr"] # jms corr
    if obj.isData:
        obj.tmp_fatjets["FatJet_msoftdrop_alt"] = obj.tmp_fatjets["FatJet_msoftdrop_altnosmear"]
        obj.tmp_fatjets["FatJet_msoftdrop"]     = obj.tmp_fatjets["FatJet_msoftdrop_nosmear"] # jms corr and no jmr
        #obj.tmp_fatjets["FatJet_msoftdrop"]     = obj.tmp_fatjets["FatJet_msoftdrop_altnosmear"] # for testing
    else:  #if not obj.isData:
        # smear mc and calc up/down jmr
        smear_corr_vals = jmr_helper(obj, opt='corr')
        smear_vals = jmr_helper(obj)
        jms_vals = {k:v for k,v in zip(['nom','down','up'],jmsVals[obj.year]) }
        # old bad sdm
        obj.tmp_fatjets["FatJet_msoftdrop_alt"] = jms_vals['nom'] * smear_vals['nom'] \
                                                  * obj.tmp_fatjets['FatJet_corr_JER'] \
                                                  * obj.tmp_fatjets["FatJet_msoftdrop"]
        def apply_jmsjmr_corr(iters_):
            # loop through all jmr and jms variations
            for k,v in iters_[0].items():
                obj.tmp_fatjets[f"FatJet_msoftdrop{v}"] = jms_vals[k] * smear_corr_vals['nom'] \
                                                           * obj.tmp_fatjets['FatJet_corr_JER'] \
                                                           * obj.tmp_fatjets["FatJet_msoftdrop_corr"]
            for k,v in iters_[1].items():
                obj.tmp_fatjets[f"FatJet_msoftdrop{v}"] = jms_vals['nom'] * smear_corr_vals[k] \
                                                          * obj.tmp_fatjets['FatJet_corr_JER'] \
                                                          * obj.tmp_fatjets["FatJet_msoftdrop_corr"]
        # dynamically handle ak8 jec variations
        if jec_sys is None or 'ak4' in jec_sys:
            apply_jmsjmr_corr([{'nom':''},{'nom':''}])
        else:
            ud  = re.search(r'((jms|jmr)(Up|Down))', jec_sys)
            if ud is not None: # jmr or jms 
                jms_iter = {'nom':''} if 'jms' not in ud.group() else {ud.group().lstrip('jms').lower():ud.group()}
                jmr_iter = {'nom':''} if 'jmr' not in ud.group() else {ud.group().lstrip('jmr').lower():ud.group()}
                apply_jmsjmr_corr([jms_iter,jmr_iter])
                obj.tmp_fatjets[f"FatJet_msoftdrop"] = obj.tmp_fatjets[f"FatJet_msoftdrop{ud.group()}"]
                del obj.tmp_fatjets[f"FatJet_msoftdrop{ud.group()}"] 
            else: # need to handle ak8jer and ak8jes variations
                apply_jmsjmr_corr([{'nom':''},{'nom':''}])
                # keep varied pt
                obj.tmp_fatjets['FatJet_pt'] = obj.fatjets['FatJet_pt']
                # re-center msoftdrop variation
                obj.tmp_fatjets['FatJet_msoftdrop'] = obj.tmp_fatjets['FatJet_msoftdrop']* (obj.fatjets['FatJet_msoftdrop'] /obj.tmp_fatjets['FatJet_msoftdrop_nom'])
                
        
def jms_helper(obj):
    getTLVm = TLorentzVectorArray.from_ptetaphim
    sjtlv = getTLVm(obj.subjets['SubJet_pt']   * (1 - obj.subjets['SubJet_rawFactor']), 
                    obj.subjets['SubJet_eta'], obj.subjets['SubJet_phi'], 
                    obj.subjets['SubJet_mass'] * (1 - obj.subjets['SubJet_rawFactor']))
    # get the raw softdrop mass from subjets
    fjsj1, fjsj2 = obj.tmp_fatjets['FatJet_subJetIdx1'], obj.tmp_fatjets['FatJet_subJetIdx2']
    groomed_mass = (sjtlv[fjsj1] + sjtlv[fjsj2]).mass
    # w-jms correction
    genCorr  = puppicorr_gen[obj.year](obj.tmp_fatjets["FatJet_pt"])
    recoCorr = aj.fromoffsets(
        obj.tmp_fatjets['FatJet_pt'].offsets,
        np.where(
            abs(obj.tmp_fatjets["FatJet_eta"].flatten()) <= 1.3,
            puppicorr_reco0eta1p3[obj.year](obj.tmp_fatjets["FatJet_pt"].flatten()),
            puppicorr_reco1p3eta2p5[obj.year](obj.tmp_fatjets["FatJet_pt"].flatten())
        )
    )
    obj.tmp_fatjets["FatJet_msoftdrop_corr"] = groomed_mass * genCorr * recoCorr
    
def jmr_helper(obj, opt=None):
    opt = '' if opt is None else '_'+opt
    getTLVm = TLorentzVectorArray.from_ptetaphim
    k_vars = ['pt','eta','phi','mass']
    gensjtlv = getTLVm(*[obj.gensubjets[f'SubGenJetAK8_{k}'] for k in k_vars])
    # delta r ak8 jet , gen jet < .4
    # delta r gen jet , gen sj jet < 0.8
    fj_gjIdx  = argmatch(obj.tmp_fatjets['FatJet_eta'],obj.tmp_fatjets['FatJet_phi'],
                         obj.genfatjets['GenJetAK8_eta'],obj.genfatjets['GenJetAK8_phi'],
                         0.4).pad(1).fillna(-1)
    gj_gsjIdx1 = argmatch(obj.genfatjets['GenJetAK8_eta'],obj.genfatjets['GenJetAK8_phi'],
                          obj.gensubjets['SubGenJetAK8_eta'],obj.gensubjets['SubGenJetAK8_phi'],
                          0.8)
    gj_gsjIdx2 = argmatch(obj.genfatjets['GenJetAK8_eta'],obj.genfatjets['GenJetAK8_phi'],
                          obj.gensubjets['SubGenJetAK8_eta'],obj.gensubjets['SubGenJetAK8_phi'],
                          0.8, m_idx=2)
    gj_mass = (gensjtlv[gj_gsjIdx1]+gensjtlv[gj_gsjIdx2]).mass    
    gj_mass = aj.fromoffsets(
        gj_gsjIdx1.offsets, 
        np.where(
            (gj_gsjIdx1.flatten() != -1) & (gj_gsjIdx2.flatten() != -1),
            gj_mass.flatten(), 
            -1
        )
    ).pad(1).fillna(-100)
    fj_gj_mass = aj.fromoffsets(
        fj_gjIdx.offsets,
        np.where(
            fj_gjIdx.flatten() != -1,
            gj_mass[fj_gjIdx].flatten(), 
            -1
        )
    )
    fj_mass = obj.tmp_fatjets[f'FatJet_msoftdrop{opt}'].pad(1).fillna(-10.)
    fj_eta = obj.tmp_fatjets['FatJet_eta'].pad(1).fillna(10.)
    fj_pt = obj.tmp_fatjets['FatJet_pt'].pad(1).fillna(-1.)
    # jmr smearer
    dM = np.where(fj_gj_mass.flatten() == -1, -1, (fj_mass - fj_gj_mass).flatten())
    smear_vals = {}
    for i, cen_or_var in enumerate(['nom','down','up']):
        # case 1 
        smearFactor = np.where((dM != -1), # gen matched 
                               1. + ((jmrVals[obj.year][i] - 1.) * dM / fj_mass.flatten()),
                               -1
        )
        # case 2 and case 3
        jm_res = np.where(
            (abs(fj_eta.flatten()) <= 1.3), 
            puppicorr_massReso_0eta1p3(fj_pt.flatten()), 
            puppicorr_massReso_1p3eta2p5(fj_pt.flatten())
        )
        rand = np.random.normal(0,jm_res)
        smearFactor = np.where(
            (dM == -1) & (fj_mass.flatten()>=0), # not gen matched and non-zero fj msoftdrop
            np.where(jmrVals[obj.year][i]>1., 1.+ rand * np.sqrt(np.power(jmrVals[obj.year][i],2) - 1.), 1),
            smearFactor)
        smearFactor = np.where(smearFactor*fj_mass.flatten() < .01, .01, smearFactor)
        jmr_smear = aj.fromoffsets(fj_mass.offsets, 
                                   smearFactor)
        smear_vals[cen_or_var] = jmr_smear[fj_mass>=0.0] # clean up anomolies
    return smear_vals
