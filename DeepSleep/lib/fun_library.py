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
import time
from functools import reduce
from numba import njit, jit, prange
import numpy as np
import pandas as pd
np.random.seed(0)
np.seterr(invalid='ignore')


######## Basic Calculations ########
def fill1e(one_d_array):
    return pd.DataFrame.from_records(one_d_array).values.flatten()
def old_fillne(n_d_array):
    return pd.DataFrame.from_records(n_d_array).values
    #
def fillne(nd_array):
    counts_ = np.array(nd_array.counts)
    rect_ , flat_ = np.full((len(nd_array),max(counts_)),np.nan), np.array(nd_array.flatten())
    @njit
    def cfillne(ja, o, c):
        rows, _ = o.shape
        c_i = 0
        for i in range(rows):
            o[i,:c[i]] = ja[c_i:c_i+c[i]]
            c_i += c[i]
        return o
    return cfillne(flat_, rect_, counts_)
    #
def sortbyscore(vars_, score_, cut_):
    ret_vars_ = []
    ind_=np.argsort(fillne(-score_[cut_]),axis=1)
    for var_ in vars_:
        temp_ = var_[:][cut_]
        ret_vars_.append(np.take_along_axis(fillne(temp_),ind_, axis=1))
    return ret_vars_
    #

def ak_crosscleaned(eta1,phi1,eta2,phi2, cut) : # cross cleaned for arbitrary length array inputs
    # input 4 awkward arrays
    # jet should be ph2,eta2
    from awkward import JaggedArray as aj
    import math
    c1_ = np.array(eta1.counts)
    c2_ = np.array(eta2.counts)
    out_ = np.ones(len(eta2.flatten()))
    args_ = eta1.flatten(), phi1.flatten(), c1_, eta2.flatten(), phi2.flatten(), c2_, cut, math.pi, out_
    @njit(parallel=False)
    def njit_cc(e1,p1,c1,e2,p2,c2,cut,pi,o):
        iter1 = 0
        iter2 = 0
        a = c1.size
        for i in range(a): 
            for ij in range(c2[i]):
                for ii in range(c1[i]):
                    deta = e1[ii+iter1] - e2[ij+iter2]
                    dphi = p1[ii+iter1] - p2[ij+iter2]
                    if (dphi > pi):
                        dphi = dphi - 2*pi
                    elif (dphi <= -pi):
                        dphi = dphi + 2*pi
                    dr = float((deta**2 + dphi**2)**(.5))
                    if dr <= cut:
                        o[ij+iter2] = 0
                        break
            iter1 += c1[i]
            iter2 += c2[i]
        return o
    out =  aj.fromcounts(eta2.counts, njit_cc(*args_)).astype(bool)
    return out

def argmatch(eta1,phi1,eta2,phi2, cut, m_idx=1): 
    # return index where 1 matches 2
    from awkward import JaggedArray as aj
    import math
    c1_ = np.array(eta1.counts).astype(int)
    c2_ = np.array(eta2.counts).astype(int)
    out_ = -1 * np.ones(len(eta1.flatten())) 
    args_ = eta1.flatten(), phi1.flatten(), c1_, eta2.flatten(), phi2.flatten(), c2_, cut, math.pi, out_, m_idx
    @njit(parallel=False)
    def njit_match(e1,p1,c1,e2,p2,c2,cut,pi,o,m_idx):
        iter1 = 0
        iter2 = 0
        a = int(c1.size)
        for i in range(a): 
            for j1 in range(c1[i]):
                m_iter = 0
                for j2 in range(c2[i]):
                    deta = float(e1[j1+iter1] - e2[j2+iter2])
                    dphi = float(p1[j1+iter1] - p2[j2+iter2])
                    if (dphi > pi):
                        dphi = dphi - 2*pi
                    elif (dphi <= -pi):
                        dphi = dphi + 2*pi
                    dr = float((deta**2 + dphi**2)**(.5))
                    if dr <= float(cut):
                        m_iter += 1
                        if m_iter == m_idx:
                            o[j1+iter1] = j2
                            break
            iter1 += c1[i]
            iter2 += c2[i]
        return o
    #return aj.fromoffsets(eta1.offsets, njit_match(*args_)).astype(int)
    return aj.fromcounts(eta1.counts, njit_match(*args_)).astype(int)

                                                

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
def invM(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2):
    m1sq = np.power(m1,2)
    m2sq = np.power(m2,2)
    pt1cosheta1_sq = np.power(pt1,2)*np.power(np.cosh(eta1),2)
    pt2cosheta2_sq = np.power(pt2,2)*np.power(np.cosh(eta2),2)
    try:
        E1pE22 = np.power((np.sqrt(pt1cosheta1_sq+m1sq)+np.sqrt(pt2cosheta2_sq+m2sq).T).T,2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = (np.sinh(eta1)*np.sinh(eta2).T).T
        p1dotp2 = (pt1*pt2.T).T*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - (pt1cosheta1_sq + pt2cosheta2_sq.T).T - 2*p1dotp2
    except (AttributeError):
        E1pE22 = np.power(np.sqrt(pt1cosheta1_sq+m1sq)+np.sqrt(pt2cosheta2_sq+m2sq),2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = np.sinh(eta1)*np.sinh(eta2)
        p1dotp2 = pt1*pt2*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - pt1cosheta1_sq - pt2cosheta2_sq - 2*p1dotp2
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
def invM_Em(pt1,eta1,phi1,E1,pt2,eta2,phi2,m2):
    m2sq = np.power(m2,2)
    pt1cosheta1_sq = np.power(pt1,2)*np.power(np.cosh(eta1),2)
    pt2cosheta2_sq = np.power(pt2,2)*np.power(np.cosh(eta2),2)
    try:
        E1pE22 = np.power((E1+np.sqrt(pt2cosheta2_sq+m2sq).T).T,2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = (np.sinh(eta1)*np.sinh(eta2).T).T
        p1dotp2 = (pt1*pt2.T).T*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - (pt1cosheta1_sq + pt2cosheta2_sq.T).T - 2*p1dotp2
    except (AttributeError):
        E1pE22 = np.power(E1+np.sqrt(pt2cosheta2_sq+m2sq),2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = np.sinh(eta1)*np.sinh(eta2)
        p1dotp2 = pt1*pt2*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - pt1cosheta1_sq - pt2cosheta2_sq - 2*p1dotp2
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
def compose(*functions):
    return reduce(lambda f, g: lambda x: f(g(x)), functions, lambda x: x)

def getZhbbBaseCuts(df_):
    import config.ana_cff as cfg
    base_cuts = (
        (df_['n_b_outZh']   >= 2)             &
        #(df_['n_lb_outZh']   >= 2)             &
        #(df_['n_b_outZh']   >= 1)             &
        (df_['n_ak4jets']   >= 5)             & 
        (df_['Zh_bbvLscore'] >= 0.8)          &
        ( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        (df_['passNotHadLep'] == 1) & # might add
        (df_['Zh_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['Zh_M']        >= 50)            &
        (df_['Zh_M']        <= 200))
    return base_cuts

def getFakebbvlCuts(df_):
    control_cuts = (
        (df_['n_b_outZh']   == 1)             &
        (df_['n_ak4jets']   >= 5)             & 
        (df_['Zh_bbvLscore'] >= 0.0)          &
        ( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        (df_['passNotHadLep'] == 1) & # might add
        (df_['Zh_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['Zh_M']        >= 50)            &
        (df_['Zh_M']        <= 200))
    return control_cuts

def getFakebbvlWeights(df_, obj):
    return (
        df_['weight']* np.sign(df_['genWeight'])
        * df_['topptWeight']
        * (df_['HEM_weight'] if obj.year == '2018' else 1.0)
        * (df_['lep_trigeffsf'])
        * df_['lep_sf']
        #* df_['dak8md_bbvl_sf']
        * df_['BTagWeight'] 
        * df_['puWeight']  
        * (df_['PrefireWeight'] if obj.year != '2018' else 1.0)
    )

def getZhbbWeight(df_, year):
    import config.ana_cff as cfg
    tot_weight = (df_['weight']* np.sign(df_['genWeight']) 
                  * (np.where(df_['weight']>300,0,1))
                  ###* (cfg.Lumi['Total']/cfg.Lumi[year])
                  * df_['topptWeight']
                  * (df_['HEM_weight']  if year == '2018' else 1.0 )
                  #* (df_['lep_trigeffsf'])
                  #* df_['lep_sf']
                  * df_['electron_sf']
                  * df_['muon_sf']
                  * df_['electron_trigeffsf']
                  * df_['muon_trigeffsf']
                  * df_['dak8md_bbvl_sf']
                  * df_['BTagWeight'] 
                  * df_['puWeight']  
                  * (df_['PrefireWeight'] if year != '2018' else 1.0))
    return tot_weight
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

def clop_pear_ci(k,n,cl=0.68, return_error=False):
    import scipy.stats
    alpha = 1-cl
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """
    lo = scipy.stats.beta.ppf(alpha/2, k, n-k+1)
    hi = scipy.stats.beta.ppf(1 - alpha/2, k+1, n-k)
    #lo, hi = map(np.nan_to_num,[lo, hi])
    if return_error: return [np.nan_to_num(abs(lo-(k/n))), np.nan_to_num(abs(hi-(k/n)), nan=1.)]
    return [np.nan_to_num(lo),np.nan_to_num(hi, nan=1.)]
    

def getLaLabel(str_):
    if '_201' in str_:
        str_ = str_.split('_201')[0] # get rid of year suffix
    la_str = ''
    col_str= ''
    la_col_map = {
        'data_obs':       [r'Data',
                           'black'],
        'ttZ':            [r't$\mathregular{\bar{t}}$Z',
                           'tab:blue'],
        'ttH':            [r't$\mathregular{\bar{t}}$H',
                           'gold'],
        'ttH_Hbb':        [r't$\mathregular{\bar{t}}$Htobb',
                           'gold'],
        'ttH_Hnonbb':     [r't$\mathregular{\bar{t}}$HtoNonbb',
                           'indianred'],
        'ttZ_Zbb':        [r't$\mathregular{\bar{t}}$Ztobb',
                           'blue'],
        'ttZ_Zqq':        [r't$\mathregular{\bar{t}}$Ztoqq',
                           'darkgreen'],
        'ttZ_Zllnunu':    [r't$\mathregular{\bar{t}}$Ztollnunu',
                           'olive'],

        'new_ttZbb':            [r't$\mathregular{\bar{t}}$Ztobb',
                                 'cyan'],

        'ttZbb':            [r't$\mathregular{\bar{t}}$Ztobb',
                             'blue'],
        'ttHbb':            [r't$\mathregular{\bar{t}}$Htobb',
                             'gold'],
        #'ttX':              [r't($\mathregular{\bar{t}}$)X',
        'ttX':              [r't$\mathregular{\bar{t}}$t$\mathregular{\bar{t}}$, t$\mathregular{\bar{t}}\gamma$, t$\mathregular{\bar{t}}$W tHW, tHq',
                             'tab:red'],
        'single_t':         [r'single top',
                             'tab:brown'],
        'TTBar':            [r't$\mathregular{\bar{t}}+\mathregular{LF}$, t$\mathregular{\bar{t}}+$c$\mathregular{\bar{c}}$',
                             'tab:orange'],
        'ttbb':        [r't$\mathregular{\bar{t}}+$b$\mathregular{\bar{b}}$',
                         'tab:green'],
        'tt_B':        [r't$\mathregular{\bar{t}}+$b$\mathregular{\bar{b}}$',
                             'tab:green'],
        'tt_bb':        [r't$\mathregular{\bar{t}}+$b$\mathregular{\bar{b}}$',
                             'tab:green'],
        'tt_2b':        [r't$\mathregular{\bar{t}}+$2b',
                             'tab:purple'],
        'VJets':            [r'V$+$jets',
                             'cyan'],
        'other':            ['VV',
                             'green'],
        'rare':            ['Rare',
                            'orange'],
        'tZq_had':         ['tZq_had',
                            'red'],
        'TWQ':             ['tZq_had',
                            'blue'],
        'THW':             ['tHW',
                            'green'],
        #
        'ttZ_genm_Zbb':    [r't$\mathregular{\bar{t}}$Z_genm_Zbb',
                            'blue'],
        'ttZ_genm_Zbb_bb':    [r't$\mathregular{\bar{t}}$Z$\rightarrow \mathrm{b\bar{b}}$ ${}_{\mathrm{GEN\;matched}}$',
                            'tab:blue'],
        'ttZ_genm_Zbb_b':    [r't$\mathregular{\bar{t}}$Z_genm_Zbb&b',
                            'purple'],
        'ttZ_genm_Zbb_nob':    [r't$\mathregular{\bar{t}}$Z_genm_Zbb&nob',
                            'tab:brown'],
        'ttH_genm_Hbb':    [r't$\mathregular{\bar{t}}$H_genm_Hbb',
                            'gold'],
        'ttH_genm_Hbb_bb':    [r't$\mathregular{\bar{t}}$H$\rightarrow \mathrm{b\bar{b}}$ ${}_{\mathrm{GEN\;matched}}$',
                           'magenta'],
        'ttH_genm_Hbb_b':    [r't$\mathregular{\bar{t}}$H_genm_Hbb&b',
                           'tab:red'],
        'ttH_genm_Hbb_nob':    [r't$\mathregular{\bar{t}}$H_genm_Hbb&nob',
                            'black'],
        'ttZ_notgenm_Zbb':    [r't$\mathregular{\bar{t}}$Z${}_\mathrm{other}$',
                            'cyan'],
        'ttH_notgenm_Hbb':    [r't$\mathregular{\bar{t}}$H${}_\mathrm{other}$',
                            'tab:red'],
        'TTZ':             [r't$\mathregular{\bar{t}}$Z', 
                            'blue'],
        'TTZ_bb':          [r't$\mathregular{\bar{t}}$Ztobb_ded',
                            'orange'],
        'TTZH':            [r't$\mathregular{\bar{t}}$Z/H',
                            'blue'],
        'TTZH_GenMatch':   [r't$\mathregular{\bar{t}}$Z/H_genMatchedZHbb',
                            'indianred'],
        'TTZH_noGenMatch': [r't$\mathregular{\bar{t}}$Z/H_nogenMatchedZHbb',
                            'black'],
        'TTZH_genZbb':     [r't$\mathregular{\bar{t}}$Z/H_genMatchedZbb',
                            'blue'],
        'TTZ_genZbb':      [r't$\mathregular{\bar{t}}$Ztobb_ded_genMatched',
                            'orange'],
        'TTZH_genZqq':     [r't$\mathregular{\bar{t}}$Z/H_genMatchedZqq',
                            'darkgreen'],
        'TTZH_genHbb':     [r't$\mathregular{\bar{t}}$Z/H_genMatchedHbb',
                            'gold'],
        'TTZH_Zbb':        [r't$\mathregular{\bar{t}}$Ztobb',
                            ':blue'],
        'TTZH_Zqq':        [r't$\mathregular{\bar{t}}$Ztoqq',
                            'darkgreen'],
        'TTZH_Hbb':        [r't$\mathregular{\bar{t}}$Htobb',
                            'gold'],
        'DY':              ['Drell-Yan',
                            ':orange'],
        'DiBoson':         ['VV',
                            ':olive'],
        'TriBoson':        ['VVV',
                            'pink'],
        'TTX':             [r't($\mathregular{\bar{t}}$)X',
                            'tab:red'],
        'TTBarLep':        [r't$\mathregular{\bar{t}}$Lep',
                            'green'],
        'TTBarSemi_pow':        [r't$\mathregular{\bar{t}}$Semi_pow',
                            'orange'],
        'TTBarDi_pow':        [r't$\mathregular{\bar{t}}$Di_pow',
                            'orange'],
        'TTBarLep_bb':     [r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$',
                            'pink'],
        'TTBarLep_pow_b':     [r't$\mathregular{\bar{t}}$+b_pow',
                            'green'],
        'TTBarLep_pow_2b':     [r't$\mathregular{\bar{t}}$+2b_pow',
                            'indianred'],
        'TTBarLep_pow_bb':     [r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$_pow',
                            'orange'],
        'TTbbHad_pow':       [r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$_dedpow',
                            'red'],
        'TTbbSemi_pow':       [r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$_dedpow',
                            'red'],
        'TTbbDi_pow':       [r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$_dedpow',
                            'red'],
        'TTBarLep_nobb':   [r't$\mathregular{\bar{t}}$',
                            'green'],
        'TTBarLep_pow_nobb':   [r't$\mathregular{\bar{t}}$_pow',
                                'brown'],
        'TTBarLep_bbGenMatch': [r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$_genMatched',
                             'grey'],
        'TTBarLep_nobbGenMatch':[r't$\mathregular{\bar{t}}$+b$\mathregular{\bar{b}}$_nogenMatched',
                             'green'],
        'TTBarHad':        [r't$\mathregular{\bar{t}}$Had',
                            'brown'],
        'TTBarHad_pow':        [r't$\mathregular{\bar{t}}$Had',
                            'brown'],
        'WJets':           [r'W$+$jets',
                            'cyan'],
        'ZJets':           [r'Z$+$jets',
                            'gray'],
        'QCD':             [r'QCD',
                            'purple']
    }
    la_str, col_str = la_col_map.get(str_,[str_,'k'])
    return [la_str, col_str]

# decorator to display the time it takes to run function
def t2Run(func):
    def wrapper(*args,**kwargs):
        start  = time.perf_counter()
        out = func(*args, **kwargs)
        finish = time.perf_counter()
        print(f'\nTime to finish {func.__name__}: {finish-start:.2f}\n')
        return out
    return wrapper

# decorator to save figures to given pdf


def import_mpl_settings(i=1):
    import matplotlib.pyplot as plt
    plt.rc("font", size=10, family="sans-serif", **{"sans-serif" : [u'TeX Gyre Heros', u'Helvetica', u'Arial']})
    plt.rc("xaxis", labellocation='right')
    plt.rc("yaxis", labellocation='top')
    plt.rc("legend", fontsize=10, scatterpoints=2, numpoints=1, borderpad=0.15, labelspacing=0.3,
           handlelength=0.7, handletextpad=0.25, handleheight=0.7, columnspacing=0.6,
           fancybox=False, edgecolor='none', borderaxespad=1)
    plt.rc("savefig", dpi=200)
    plt.rc("figure", figsize=(3.375*i, 3.375*(6./8.)*i), dpi=200)
    plt.rc("text.latex", preamble='\n'.join([r"\usepackage{amsmath}",
                                             r"\usepackage{helvet}",
                                             r"\usepackage{sansmath}",
                                             r"\sansmath"]))
    

def upperlefttext(s):
    trans = gca().transAxes + matplotlib.transforms.ScaledTranslation(3/72, -3/72, gcf().dpi_scale_trans)
    return text(0, 1, s, transform=trans, ha='left', va='top')

def CMSlabel(fig=None, ax=None, opt=None, altax=None, lumi=None):
    import matplotlib.pyplot as plt
    from matplotlib import transforms
    if fig is None:
        fig = plt.gcf()
    if ax is None:
        ax = plt.gca()
    if altax is None:
        ax0 = ax1 = ax
    else:
        ax0 = ax
        ax1 = altax
    if opt is None:
        opt = 'Preliminary'
    if lumi is None:
        lumi = 137
    lumi = f'{lumi:.1f}' if float(lumi) < 100 else str(lumi)
    trans = ax0.transAxes + transforms.ScaledTranslation(0/72, 3/72, fig.dpi_scale_trans)
    ax0.text(0, 1, rf'\textbf{{CMS}} {{\footnotesize \textit{{{opt}}}}}', usetex=True,
             transform=trans, ha='left', va='baseline')
    trans = ax1.transAxes + transforms.ScaledTranslation(0/72, 3/72, fig.dpi_scale_trans)
    ax1.text(1, 1, rf"{{\footnotesize ${{{lumi}}}\,\mathrm{{fb}}^{{\text{{-1}}}}$ (13 TeV)}}", usetex=True,
             transform=trans, ha='right', va='baseline')

def make_error_boxes(ax, xdata, ydata, xerror, yerror,  facecolor='r',
                     edgecolor='None', alpha=0.0, hatch=10*'X', label=''):
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    errorboxes = []
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe, y - ye), 2*xe, 2*ye) if type(ye) is not np.ndarray else Rectangle((x - xe, y - ye[0]), 2*xe, ye[0]+ye[1])
        errorboxes.append(rect)
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor, hatch=hatch, label=label, zorder=1.5)
    # Add collection to axes
    ax.add_collection(pc)

def save_pdf(pdf_name = 'dummy.pdf'):
    import matplotlib.backends.backend_pdf as matpdf 
    import matplotlib.pyplot as plt
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)',
        shell=True).decode().strip('\n')+'DeepSleep/')
    if pdf_name == '':
        raise NameError("Yo, dummy, put a filename here ''")
    def inner (func):
        def wrapper(*args,**kwargs):
            pdf = matpdf.PdfPages(f"{sys.path[1]}pdf/{pdf_name}")  
            func(*args, **kwargs) # doesnt return anything
            print(f"Saving figures to: {sys.path[1]}pdf/{pdf_name}")
            for fig_ in range(1, plt.gcf().number+1):
                pdf.savefig( fig_ )
            pdf.close()
            plt.close('all')

        return wrapper
    return inner
