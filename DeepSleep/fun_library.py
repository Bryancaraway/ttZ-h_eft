#                      #  
##                    ##
########################                               
### Library for      ###
### TTZ/H utility    ###
### functions        ###
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
np.random.seed(0)


######## Basic Calculations ########
def fill1e(one_d_array):
    return pd.DataFrame.from_records(one_d_array).values.flatten()
def fillne(n_d_array):
    return pd.DataFrame.from_records(n_d_array).values
    #
def sortbyscore(vars_, score_, cut_):
    ret_vars_ = []
    ind_=np.argsort(fillne(-score_[cut_]),axis=1)
    for var_ in vars_:
        temp_ = var_[:][cut_]
        ret_vars_.append(np.take_along_axis(fillne(temp_),ind_, axis=1))
    return ret_vars_
    #
def deltaR(eta1,phi1,eta2,phi2):
    try:
        deta = np.subtract(eta1,eta2.T).T
        dphi = np.subtract(phi1,phi2.T).T
    except (AttributeError) :
        deta = eta1 - eta2
        dphi = phi1 - phi2
    #
    dphi[((dphi > math.pi)   & (dphi != np.nan))] = dphi[((dphi > math.pi)   & (dphi != np.nan))] - 2*math.pi
    dphi[((dphi <= -math.pi) & (dphi != np.nan))] = dphi[((dphi <= -math.pi) & (dphi != np.nan))] + 2*math.pi
    #
    delta_r = np.sqrt(np.add(np.power(deta,2),np.power(dphi,2)))
    return delta_r
    #
def invM(pt1,eta1,phi1,pt2,eta2,phi2):
    try:
        pt1pt2 = (pt1*pt2.T).T
        cosheta1eta2 = np.cosh((eta1-eta2.T).T)
        cosphi1phi2  = np.cos(deltaPhi(phi1,phi2))
    except (AttributeError) :
        pt1pt2 = pt1*pt2
        cosheta1eta2 = np.cosh(eta1-eta2)
        cosphi1phi2  = np.cos(deltaPhi(phi1,phi2))
    #
    invm2 = 2*pt1pt2*(cosheta1eta2-cosphi1phi2)
    return np.sqrt(invm2)
def invM_sdM(pt1,eta1,phi1,m1,pt2,eta2,phi2,E2):
    m1sq = np.power(m1,2)
    pt1cosheta1_sq = np.power(pt1,2)*np.power(np.cosh(eta1),2)
    pt2cosheta2_sq = np.power(pt2,2)*np.power(np.cosh(eta2),2)
    try:
        E1pE22 = np.power((np.sqrt(pt1cosheta1_sq+m1sq)+E2.T).T,2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = (np.sinh(eta1)*np.sinh(eta2).T).T
        p1dotp2 = (pt1*pt2.T).T*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - (pt1cosheta1_sq + pt2cosheta2_sq.T).T - 2*p1dotp2
    except (AttributeError):
        E1pE22 = np.power(np.sqrt(pt1cosheta1_sq+m1sq)+E2,2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = np.sinh(eta1)*np.sinh(eta2)
        p1dotp2 = pt1*pt2*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - pt1cosheta1_sq - pt2cosheta2_sq - 2*p1dotp2
    return np.sqrt(invm2)
def invM_E(pt1,eta1,phi1,E1,pt2,eta2,phi2,E2):
    pt1cosheta1_sq = np.power(pt1,2)*np.power(np.cosh(eta1),2)
    pt2cosheta2_sq = np.power(pt2,2)*np.power(np.cosh(eta2),2)
    cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
    try:
        sinheta1Xsinheta2 = (np.sinh(eta1)*(np.sinh(eta2)).T).T
        p1dotp2 = (pt1*pt2.T).T*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = np.power(E1+E2.T,2).T - (pt1cosheta1_sq + (pt2cosheta2_sq + 2*p1dotp2).T).T
    except (AttributeError) :
        sinheta1Xsinheta2 = np.sinh(eta1)*np.sinh(eta2)
        p1dotp2 = pt1*pt2*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = np.power(E1+E2,2) - pt1cosheta1_sq - pt2cosheta2_sq - 2*p1dotp2 
    return np.sqrt(invm2)
    #
def deltaPhi(phi1, phi2):
    try:
        dphi = np.subtract(phi1,phi2.T).T
    except (AttributeError) :
        dphi = phi1-phi2
    dphi[((dphi > math.pi)   & (dphi != np.nan))] = dphi[((dphi > math.pi)   & (dphi != np.nan))] - 2*math.pi
    dphi[((dphi <= -math.pi) & (dphi != np.nan))] = dphi[((dphi <= -math.pi) & (dphi != np.nan))] + 2*math.pi
    
    return dphi

def calc_mtb(pt_, phi_, m_pt, m_phi):
    try:
        mtb2 =  2*(m_pt*pt_.T).T*(1 - np.cos(deltaPhi(m_phi,phi_)))
    except (AttributeError):
         mtb2 =  2*m_pt*pt_*(1 - np.cos(deltaPhi(phi_, m_phi)))
    return np.sqrt(mtb2)
    #
def calc_SandA(pt_,eta_,phi_): # sphericity and aplanarity
    S_      = np.zeros((pt_.shape[0],3,3))
    pxyz_   = np.array([pt_*np.cos(phi_),pt_*np.sin(phi_),pt_*np.sinh(eta_)])
    p2_sum  = np.nansum(np.power(pt_*np.cosh(eta_),2),axis=1)
    for i_ in range(0,3):
        for j_ in range(0,3):
            S_[:,i_,j_] = np.nansum(pxyz_[i_]*pxyz_[j_], axis = 1)/p2_sum
    #
    eig_= -np.sort(-(np.linalg.eig(S_[:])[0]), axis = 1) # calc eig vals and sort by descending
    s_ = (3/2)*(eig_[:,1]+eig_[:,2])
    a_ = (3/2)*eig_[:,2]
    return s_, a_
########### Auxihilary Functions #############
def getZHbbBaseCuts(dict_):
    base_cuts = (
        (dict_['ak8']['n_nonHbb'] >= 2)    &
        (dict_['ak8']['nhbbFatJets'] > 0)  &
        (dict_['ak8']['H_pt']       >= 200)&
        (dict_['ak8']['H_M']         > 50) &
        (dict_['ak8']['H_M']         < 200))
    return base_cuts
#
def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)
#
def StackedHisto(df_, kinem_, range_, xlabel_, n_bins=20):
    from matplotlib import rc
    #
    fontsize = 12    
    rc("savefig", dpi=250)
    #rc("figure", figsize=(3.375, 3.375*(6./8.)), dpi=250)                                                            
    #rc("text.latex", preamble=r"\usepackage{amsmath}")                                                             
    #rc('font',**{'family':'serif','serif':['Times']})
    #rc("hatch", linewidth=0.0) 
    #
    h        = []
    h_sum    = []
    w        = []
    integral = []
    labels   = []
    colors   = []
    #

    bins = np.arange(range_[0],range_[-1]+((range_[-1])/n_bins) , (range_[-1])/n_bins)
    #                                                                                                
    for i_, key_ in enumerate(df_.keys()):
        if (kinem_ in df_[key_]['df'].keys()):
            kinem     = df_[key_]['df'][kinem_]
        elif (kinem_ in df_[key_]['val'].keys()):
            kinem     = df_[key_]['val'][kinem_]
        elif (kinem_ in df_[key_]['valRC'].keys()):
            kinem     = df_[key_]['valRC'][kinem_]
        try:
            if (kinem_ in df_[key_]['ak8'].keys()):
                kinem = df_[key_]['ak8'][kinem_]
        except:
            pass
        base_cuts = (
            #(df_[key_]['val']['NN'] <  .85) & 
            (df_[key_]['val']['NN'] >= 0.0) & 
            
            (df_[key_]['ak8']['H_pt'] >= 200)       &
            #(df_[key_]['ak8']['H_pt'] < 300)       &
            #(df_[key_]['val']['genZHpt'] >= 300)       &
            #(df_[key_]['val']['genZHpt'] < 350)       &
            #(df_[key_]['ak8']['n_H_sj_btag'] == 2) &
            #(df_[key_]['ak8']['n_b_Hbb'] == 2) &
            
            (df_[key_]['ak8']['n_nonHbb'] >= 2)    &
            (df_[key_]['ak8']['nhbbFatJets'] > 0)  &
            (df_[key_]['ak8']['H_M']         > 50) &  
            (df_[key_]['ak8']['H_M']         < 200)& 
            (df_[key_]['val']['MET_pt']      >= 0))# &
        if ('_GenMatch' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == True)
        if ('_noGenMatch' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == False) & (df_[key_]['val']['matchedGen_Zqq'] == False) 
        if ('_genZbb' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['matchedGen_Zbb'] == True)
        if ('_genZqq' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['matchedGen_Zqq'] == True)
        if ('_genHbb' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['matchedGen_Hbb'] == True)
        if ('_Zbb' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['Zbb'] == True)
        if ('_Zqq' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['Zqq'] == True)
        if ('_Hbb' in key_):
            base_cuts = base_cuts & (df_[key_]['val']['Hbb'] == True)
        ########
        h.append( np.clip(kinem[base_cuts], bins[0], bins[-1]))
        w.append( df_[key_]['val']['weight'][base_cuts] * np.sign(df_[key_]['val']['genWeight'][base_cuts]) * (137/41.9)*\
                  df_[key_]['val']['BTagWeight'][base_cuts]*df_[key_]['val']['puWeight'][base_cuts]*df_[key_]['val']['PrefireWeight'][base_cuts])
        n_, bins_,_ = plt.hist(h[i_], weights=w[i_])
        integral.append( sum(n_[:]))
        la_label, color = kFit.getLaLabel(key_)
        labels.append( la_label + ' ({0:3.1f})'.format(integral[i_]))
        colors.append( color)
        plt.close('all')
    #          
    fig, ax = plt.subplots()
    fig.subplots_adjust(
        top=0.88,
        bottom=0.11,
        left=0.11,
        right=0.88,
        hspace=0.2,
        wspace=0.2
    )
    #fig.subplots_adjust(left=0.03, right=0.97, bottom=0.05, top=0.92)
    perc_plot = False
    n_, bins_, patches_ = ax.hist(h,
                                  bins=bins, stacked=False,# fill=True,
                                  #range=range_,
                                  histtype='step',
                                  density=False,
                                  #linewidth=0,
                                  weights= w if not perc_plot else np.divide(w,integral),
                                  color  = colors,
                                  label= labels)
    #
    #ax.grid(True)
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    fig.text(0.12,0.89, r"$\bf{CMS}$ $Simulation$",   fontsize = fontsize)
    fig.text(0.64,0.89, r'137 fb$^{-1}$ (13 TeV)', fontsize = fontsize)
    plt.xlabel(xlabel_, fontsize = fontsize)
    plt.ylabel(f"{'%' if perc_plot else 'Events'} / {range_[-1]/n_bins:.2f}", 
               fontsize = fontsize)
    plt.xlim(range_)
    if not perc_plot: plt.yscale('log')
    plt.grid(True)
    #plt.setp(patches_, linewidth=0)
    plt.legend(framealpha = 1, ncol=2)
    plt.savefig('money_pdf/moneyplot'+xlabel_+'_.pdf', dpi = 300)
    plt.show()
    plt.close(fig)
    #    
# 
def calc_Kappa(df_,nn_range):
    sig = 'TTZH'
    bkg = 'TTBarLep'
    suffix = '_2017'
    def bin_counts(key_):
        base_cuts = (
            (df_[key_]['ak8']['n_nonHbb'] >= 2)    &
            (df_[key_]['ak8']['nhbbFatJets'] > 0)  &
            (df_[key_]['ak8']['H_M']         > 50) &
            (df_[key_]['ak8']['H_M']         < 200))
        if ('TTZH' in key_ ):
            base_cuts = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == True)
        w = df_[key_]['val']['weight'][base_cuts] * np.sign(df_[key_]['val']['genWeight'][base_cuts]) * (137/41.9)
        a_cut = (
            (df_[key_]['ak8']['H_M'] > 70) & (df_[key_]['ak8']['H_M'] < 145) &
            (df_[key_]['val']['NN'] >= .95))
        az_cut, ah_cut = a_cut & (df_[key_]['ak8']['H_M'] <= 100), a_cut & (df_[key_]['ak8']['H_M'] > 100)
        b_cut = (
            (((df_[key_]['ak8']['H_M']> 50) & (df_[key_]['ak8']['H_M'] <= 70)) | ((df_[key_]['ak8']['H_M'] >= 145) & (df_[key_]['ak8']['H_M'] < 200))) &
            (df_[key_]['val']['NN'] >= .95))
        bl_cut, br_cut = b_cut & (df_[key_]['ak8']['H_M'] <= 70), b_cut & (df_[key_]['ak8']['H_M'] >= 145)
        c_cut = (
            (df_[key_]['ak8']['H_M'] > 70) & (df_[key_]['ak8']['H_M'] < 145) &
            (df_[key_]['val']['NN'] < nn_range[1]) & (df_[key_]['val']['NN'] >= nn_range[0]))
        cz_cut, ch_cut = c_cut & (df_[key_]['ak8']['H_M'] <= 100), c_cut & (df_[key_]['ak8']['H_M'] > 100)
        d_cut = (
            (((df_[key_]['ak8']['H_M']> 50) & (df_[key_]['ak8']['H_M'] <= 70)) | ((df_[key_]['ak8']['H_M'] >= 145) & (df_[key_]['ak8']['H_M'] < 200))) &
            (df_[key_]['val']['NN'] < nn_range[1]) & (df_[key_]['val']['NN'] >= nn_range[0]))
        dl_cut, dr_cut = d_cut & (df_[key_]['ak8']['H_M'] <= 70), d_cut & (df_[key_]['ak8']['H_M'] >= 145)
        return [w[a_cut].sum(), w[az_cut].sum(), w[ah_cut].sum(),
                w[b_cut].sum(), w[bl_cut].sum(), w[br_cut].sum(),
                w[c_cut].sum(), w[cz_cut].sum(), w[ch_cut].sum(),
                w[d_cut].sum(), w[dl_cut].sum(), w[dr_cut].sum()]
    #
    a_sig, az_sig, ah_sig, b_sig, bl_sig, br_sig, c_sig, cz_sig, ch_sig, d_sig, dl_sig, dr_sig = bin_counts(sig+suffix)
    a_sig_err, az_sig_err, ah_sig_err, b_sig_err, bl_sig_err, br_sig_err, c_sig_err, cz_sig_err, ch_sig_err, d_sig_err, dl_sig_err, dr_sig_err = np.sqrt([
        a_sig, az_sig, ah_sig, b_sig, bl_sig, br_sig, c_sig, cz_sig, ch_sig, d_sig, dl_sig, dr_sig])
    a_bkg, az_bkg, ah_bkg, b_bkg, bl_bkg, br_bkg, c_bkg, cz_bkg, ch_bkg, d_bkg, dl_bkg, dr_bkg = bin_counts(bkg+suffix)
    a_bkg_err, az_bkg_err, ah_bkg_err, b_bkg_err, bl_bkg_err, br_bkg_err, c_bkg_err, cz_bkg_err, ch_bkg_err, d_bkg_err, dl_bkg_err, dr_bkg_err = np.sqrt([
        a_bkg, az_bkg, ah_bkg, b_bkg, bl_bkg, br_bkg, c_bkg, cz_bkg, ch_bkg, d_bkg, dl_bkg, dr_bkg])
    #
    k_bkg     = (b_bkg * c_bkg)/(a_bkg * d_bkg)
    k_bkg_err = abs(k_bkg)*np.sqrt(np.power(a_bkg_err/a_bkg,2) + np.power(b_bkg_err/b_bkg,2) + np.power(c_bkg_err/c_bkg,2) + np.power(d_bkg_err/d_bkg,2)) 
    a_eff     = (a_sig/a_bkg)
    a_eff_err = abs(a_eff)*np.sqrt(np.power(a_sig_err/a_sig,2) + np.power(a_bkg_err/a_bkg,2))
    c_eff = (c_sig/c_bkg)
    c_eff_err = abs(c_eff)*np.sqrt(np.power(c_sig_err/c_sig,2) + np.power(c_bkg_err/c_bkg,2))
    print('\nFor NN score range: {0:1.2f}-->{1:1.2f}\n'\
          '=========================================\n'\
          'Bin_A sig/bkg efficiency: {2:2.2f} +/- {3:2.2f}%\n'\
          'Bin_C sig/bkg efficiency: {4:2.2f} +/- {5:2.2f}%\n'\
          'Kappa: {6:1.2f} +/- {7:1.2f}\n'.format(nn_range[0], nn_range[1], a_eff*100, a_eff_err*100, c_eff*100, c_eff_err*100, k_bkg, k_bkg_err))
    kz_bkg = (bl_bkg * br_bkg  * cz_bkg)/(az_bkg * dl_bkg * dr_bkg)
    kz_bkg_err = abs(kz_bkg)*np.sqrt(np.power(az_bkg_err/az_bkg,2) + np.power(bl_bkg_err/bl_bkg,2) + np.power(br_bkg_err/br_bkg,2) + 
                                     np.power(cz_bkg_err/cz_bkg,2) + np.power(dl_bkg_err/dl_bkg,2) + np.power(dr_bkg_err/dr_bkg,2))
    kh_bkg = (bl_bkg * br_bkg  * ch_bkg)/(ah_bkg * dl_bkg * dr_bkg)
    kh_bkg_err = abs(kh_bkg)*np.sqrt(np.power(ah_bkg_err/ah_bkg,2) + np.power(bl_bkg_err/bl_bkg,2) + np.power(br_bkg_err/br_bkg,2) + 
                                     np.power(ch_bkg_err/ch_bkg,2) + np.power(dl_bkg_err/dl_bkg,2) + np.power(dr_bkg_err/dr_bkg,2))
    az_eff, ah_eff = (az_sig/az_bkg), (ah_sig/ah_bkg)
    az_eff_err, ah_eff_err = abs(az_eff)*np.sqrt(np.power(az_sig_err/az_sig,2) + np.power(az_bkg_err/az_bkg,2)), abs(ah_eff)*np.sqrt(np.power(ah_sig_err/ah_sig,2) + np.power(ah_bkg_err/ah_bkg,2))
    cz_eff, ch_eff = (cz_sig/cz_bkg), (ch_sig/ch_bkg)
    cz_eff_err, ch_eff_err = abs(cz_eff)*np.sqrt(np.power(cz_sig_err/cz_sig,2) + np.power(cz_bkg_err/cz_bkg,2)), abs(ch_eff)*np.sqrt(np.power(ch_sig_err/ch_sig,2) + np.power(ch_bkg_err/ch_bkg,2))
    print(
        '=========================================\n'\
        'Bin_AZ sig/bkg efficiency: {2:2.2f} +/- {3:2.2f}%\n'\
        'Bin_CZ sig/bkg efficiency: {4:2.2f} +/- {5:2.2f}%\n'\
        'KappaZ: {6:1.2f} +/- {7:1.2f}\n'\
        '--\t--\t--\t--\t--\t--\n'\
        'Bin_AH sig/bkg efficiency: {8:2.2f} +/- {9:2.2f}%\n'\
        'Bin_CH sig/bkg efficiency: {10:2.2f} +/- {11:2.2f}%\n'\
        'KappaH: {12:1.2f} +/- {13:1.2f}\n'.format(
            nn_range[0], nn_range[1],
            az_eff*100, az_eff_err*100, 
            cz_eff*100, cz_eff_err*100, 
            kz_bkg, kz_bkg_err, 
            ah_eff*100, ah_eff_err*100, 
            ch_eff*100, ch_eff_err*100, 
            kh_bkg, kh_bkg_err))
    #
    #return k_bkg, k_bkg_err, c_eff, c_eff_err
    bdl_bkg , bdl_bkg_err = bl_bkg/dl_bkg, abs(bl_bkg/dl_bkg)*np.sqrt(np.power(bl_bkg_err/bl_bkg,2) + np.power(dl_bkg_err/dl_bkg,2))
    bdr_bkg , bdr_bkg_err = br_bkg/dr_bkg, abs(br_bkg/dr_bkg)*np.sqrt(np.power(br_bkg_err/br_bkg,2) + np.power(dr_bkg_err/dr_bkg,2))
    acz_bkg , acz_bkg_err = az_bkg/cz_bkg, abs(az_bkg/cz_bkg)*np.sqrt(np.power(az_bkg_err/az_bkg,2) + np.power(cz_bkg_err/cz_bkg,2))
    ach_bkg , ach_bkg_err = ah_bkg/ch_bkg, abs(ah_bkg/ch_bkg)*np.sqrt(np.power(ah_bkg_err/ah_bkg,2) + np.power(ch_bkg_err/ch_bkg,2))
    #
    bl_eff, br_eff = (bl_sig/bl_bkg), (br_sig/br_bkg)
    bl_eff_err, br_eff_err = abs(bl_eff)*np.sqrt(np.power(bl_sig_err/bl_sig,2) + np.power(bl_bkg_err/bl_bkg,2)), abs(br_eff)*np.sqrt(np.power(br_sig_err/br_sig,2) + np.power(br_bkg_err/br_bkg,2))
    #
    return [
        [bdl_bkg,acz_bkg,ach_bkg,bdr_bkg],
        [bdl_bkg_err,acz_bkg_err,ach_bkg_err,bdr_bkg_err],
        [bl_eff,az_eff,ah_eff,br_eff],
        [bl_eff_err,az_eff_err,ah_eff_err,br_eff_err]
    ]
