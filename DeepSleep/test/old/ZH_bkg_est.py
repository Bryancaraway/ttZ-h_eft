#                      #  
##                    ##
########################                               
### TTZ/H, Z/H to bb ###
### signal and bkg   ###                               
### SR estimation    ###                               
########################                               
### written by:      ###                               
### Bryan Caraway    ###                               
########################                               
##                    ##                                 
#                      #

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
np.random.seed(1)
##


def Bkg_Est(files_, samples_, outDir_, overlap_ = cfg.ZHbbFitoverlap, n_bins = 14):
    from fun_library import calc_Kappa
    df = kFit.retrieveData(files_, samples_, outDir_, getgen_=False, getak8_=True)
    sig = 'TTZH'
    bkg = 'TTBarLep'
    rare = ['TriBoson','DiBoson','QCD','WJets','DY']
    ttX  = ['TTX','ttZqq']
    suf = '_2017'
    NN_cut = .96
    #
    m_bins = np.arange(50.0,200.0,(200.0-50.0)/n_bins)
    #
    def applyCuts(df_,key_):
        base_cuts =(
            (df_[key_]['ak8']['n_nonHbb'] >= 2)    &
            (df_[key_]['ak8']['nhbbFatJets'] > 0)  &
            (df_[key_]['ak8']['H_M']         > 50) &
            (df_[key_]['ak8']['H_M']         < 200))
        if (sig in key_): # look at the sig/bkg ratios only in the signal region
            base_cuts = base_cuts & (df_[key_]['val']['NN'] >= NN_cut)
            _genm   = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == True)
            _nogenm = base_cuts & (df_[key_]['val']['matchedGen_ZHbb'] == False)
            _genZ   = base_cuts & (df_[key_]['val']['matchedGen_Zbb'] == True)
            _genH   = base_cuts & (df_[key_]['val']['matchedGen_Hbb'] == True)
            _genZbb = base_cuts & (df_[key_]['val']['Zbb'] == True)
            _genHbb = base_cuts & (df_[key_]['val']['Hbb'] == True)
            _genZqq = base_cuts & (df_[key_]['val']['Zqq'] == True)
            sig_genm, sig_genmw     = df_[key_]['ak8']['H_M'][_genm],   (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genm  ]
            sig_nogenm, sig_nogenmw = df_[key_]['ak8']['H_M'][_nogenm], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_nogenm]
            sig_genZ, sig_genZw     = df_[key_]['ak8']['H_M'][_genZ],   (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genZ  ]
            sig_genH, sig_genHw     = df_[key_]['ak8']['H_M'][_genH],   (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genH  ]
            sig_genZbb, sig_genZbbw = df_[key_]['ak8']['H_M'][_genZbb], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genZbb]
            sig_genHbb, sig_genHbbw = df_[key_]['ak8']['H_M'][_genHbb], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genHbb]
            sig_genZqq, sig_genZqqw = df_[key_]['ak8']['H_M'][_genZqq], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_genZqq]

            return [sig_genm, sig_nogenm, sig_genZ, sig_genH,
                    sig_genZbb, sig_genHbb, sig_genZqq,
                    sig_genmw, sig_nogenmw, sig_genZw, sig_genHw,
                    sig_genZbbw, sig_genHbbw, sig_genZqqw]
        else:
            _sr = base_cuts & (df_[key_]['val']['NN'] >= NN_cut)
            bkg_sr, bkg_srw, bkg_sr_gw = df_[key_]['ak8']['H_M'][_sr], (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[_sr  ], np.sign(df_[key_]['val']['genWeight'][_sr])  
            return bkg_sr, bkg_srw, bkg_sr_gw
            #
    def shapePlot(y_,y_err_,label_,ylabel_, dohist=False):
        bin_c = (m_bins[1:] + m_bins[:-1])/2
        bin_w = m_bins[1:] - m_bins[:-1]
        if dohist:
            plt.hist(x=np.repeat(bin_c[np.newaxis,:],y_.shape[0], axis=0).T, weights = y_.T,
                     bins=m_bins,
                     stacked=True,
                     histtype='step',label=label_)
        else:
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
    bkg_sr, bkg_srW, bkg_sr_gw = applyCuts(df,bkg+suf)
    bkg_df = pd.DataFrame(columns=['bkg_sr','bkg_srw','bkg'])
    sig_genm, sig_nogenm, sig_genZ, sig_genH, sig_genZbb, sig_genHbb, sig_genZqq, sig_genmw, sig_nogenmw, sig_genZw, sig_genHw, sig_genZbbw, sig_genHbbw, sig_genZqqw = applyCuts(df,sig+suf)
    for key_ in df.keys():
        if sig in key_: continue
        _sr, _srw, sr_gw = applyCuts(df,key_)
        name = key_.split('_201')[0]
        if _sr.size == 0 : continue            
        #if name == 'QCD':  continue
        if name in rare:
            name = 'Rare'
        elif name in ttX:
            name = 'ttX'
        bkg_df = bkg_df.append(pd.DataFrame({'bkg_sr':_sr,'bkg_srw':_srw, 'bkg': name}), ignore_index=True) 
    #
    bkg_df = bkg_df.append(pd.DataFrame({'bkg_sr':sig_genZqq,'bkg_srw':sig_genZqqw, 'bkg':'ttX'}), ignore_index=True)
    #print(bkg_df)
    shapePlot(*histvals(bkg_df['bkg_sr'],  val_w=bkg_df['bkg_srw'],   norm=False),    'bkg',   '')    
    shapePlot(*histvals(sig_genHbb,          val_w=sig_genHbbw,        norm=False),    'ttH',   '')    
    shapePlot(*histvals(sig_genZbb,          val_w=sig_genZbbw,        norm=False),    'ttZ',   '')    
    plt.legend()
    #plt.show()
    plt.clf()
    #
    # create data_card for chi^2 test
    # format : sample name , hist values, hist value error
    data_card = pd.DataFrame(columns=['vals','valsE'])
    _vals, _valsE = histvals(sig_genZbb, val_w=sig_genZbbw)
    data_card = data_card.append(pd.DataFrame({'vals': [_vals], 'valsE': [_valsE]}, index=['ttZbb']))
    _vals, _valsE = histvals(sig_genHbb, val_w=sig_genHbbw)
    data_card = data_card.append(pd.DataFrame({'vals': [_vals], 'valsE': [_valsE]}, index=['ttHbb']))
    for bkg_ in set(bkg_df['bkg']):
        temp_ = bkg_df[bkg_df['bkg'] == bkg_]
        _vals, _valsE = histvals(temp_['bkg_sr'], temp_['bkg_srw'])
        data_card = data_card.append(pd.DataFrame({'vals': [_vals], 'valsE': [_valsE]}, index=[bkg_]))
    data_card['rate'] = 1.0
    # calc sub values for chi2 calculation
    B_ij=np.stack(data_card['vals'].values)
    bkg_ij = np.stack(data_card['vals'].loc[['Rare', 'TTBarLep','ttX']].values)
    sig_ij = np.stack(data_card['vals'].loc[['ttZbb','ttHbb']].values)
    r_j =data_card['rate'].values
    #print(np.random.poisson(D_i,(10,D_i.size)))
    n_pdD = 35
    
    #D_i = np.average(D_i,axis=0)
    #print(D_i)
    def GLS(data_,x_):
        #
        X_     = np.matrix(x_)
        Omega_ = np.matrix(np.diag(data_))
        Y_     = np.matrix(data_).T
        #
        Cov_   = (X_ * Omega_.I * X_.T).I # transpose is flipped because X is inheirently flipped in shape
        beta_  = Cov_ * (X_ * Omega_.I * Y_)
        Corr_  = Cov_.A / np.sqrt(np.outer(np.diag(Cov_),np.diag(Cov_)))
        return beta_.A.flatten(),Cov_,Corr_
    #

    x   = []
    cov = []
    cor = []
    for _ in range(10000):
        D_i = np.random.poisson(np.sum(B_ij,axis=0)) #size=(n_pdD,len(np.sum(B_ij,axis=0))))
        x_ , cov_ , cor_  = GLS(
            D_i-np.sum(bkg_ij,axis=0),
            #B_ij)
            sig_ij)
        if (abs(x_) > 5).any(): continue
        x.append(x_)
        #
    x=np.array(x)
    print(x.shape)
    print(np.mean(x,axis=0))
    print(np.median(x,axis=0))
    print(np.std(x,axis=0).shape)
    plt.hist(x, histtype='step', label= ['ttZ','ttH'])
    plt.xlabel('$\mu$')
    plt.legend()
    plt.show()
    plt.clf()
    plt.errorbar(x=[350,350][0],y=np.median(x,axis=0)[0],xerr=[150,150][0],yerr=np.std(x,axis=0)[0], label=['ttZ','ttH'][0],
                 fmt='o',barsabove=True,capsize=5)
    plt.errorbar(x=[350,350][1],y=np.median(x,axis=0)[1],xerr=[150,150][1],yerr=np.std(x,axis=0)[1], label=['ttZ','ttH'][1],
                 fmt='o',barsabove=True,capsize=5)
    plt.legend()
    plt.ylabel('$\mu$')
    plt.xlabel('H/Z pt (GeV)')
    plt.show()
    plt.clf()
    exit()
    r_j[:2] = x
    x = r_j
    print('\n',cov)
    print('\n',cor)
    ###print(np.sum(B_ij,axis=0))
    #### assemble chi^2
    def chisqfun(_r):
        sum_r_j_B_ij = np.sum((_r*B_ij.T).T,axis=0)
        chi2 = np.nansum(np.power(D_i - sum_r_j_B_ij,2)/sum_r_j_B_ij)
        return chi2
    #### minimize chi2 function and find new rates (will need to add tolorance to rates later)
    #####
    ###import scipy.optimize as opt
    ###res = opt.minimize(chisqfun,r_j, method='BFGS')
    ####res = opt.least_squares(chisqfun,r_j)
    ###x = res.x
    ####res = opt.leastsq(chisqfun,r_j)
    #####
    ####
    print(np.sum((x*B_ij.T).T,axis=0))
    ###print('Success?',res.success)
    print(chisqfun(x))
    data_card['new_rates'] = x
    ####print(pd.DataFrame(res.jac,columns=data_card.index))
    print(data_card['new_rates'])
    print('\nEstimated Covariance Matrix')
    #print(pd.DataFrame(cov,columns=data_card.index, index=data_card.index))
    ####print(pd.DataFrame(res.hess_inv, columns=data_card.index, index=data_card.index))
    print('\nCorrelation Matrix')
    #print(pd.DataFrame(cor,columns=data_card.index, index=data_card.index))
    #print(pd.DataFrame(np.corrcoef(res.hess_inv), columns=data_card.index, index=data_card.index))
    #print(pd.DataFrame(np.cov(B_ij), columns=data_card.index, index=data_card.index))
    #
    shapePlot(D_i,np.sqrt(D_i), 'psData', '')
    #for y_,y_err_,label_ in zip((x*B_ij.T).T,np.sqrt((x*B_ij.T).T),data_card.index):
    #    print(y_.shape)
    #    print(m_bins.shape)
    #    shapePlot(y_,  y_err_,  label_,'',dohist=True)
    shapePlot((x*B_ij.T).T,  np.sqrt((x*B_ij.T).T),  data_card.index.values,'',dohist=True)
    
    plt.legend()
    plt.show()
    plt.clf()
def measSigStr(files_, samples_, outDir_):
    df_ = kFit.retrieveData(files_, samples_, outDir_, getak8_ = True)
    sig = 'TTZH'
    suf = '_2017'
    # Organize Sig and Bkg DataFrame
    sig_df = pd.DataFrame()
    bkg_df = pd.DataFrame()
    for key_ in df_.keys():
        def get_weight(cut):
            if len(cut) > 0:
                weight = (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[cut]
            else:
                weight = (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))
            return weight
        def add_toDF(df, cut=[], name=key_.split('_201')[0]):
            if len(cut) > 0:
                df = df.append(pd.DataFrame({'HZ_pt':df_[key_]['ak8']['H_pt'][cut], 'NN':df_[key_]['val']['NN'][cut], 'Weight':get_weight(cut), 'Name':name}), ignore_index=True)
            else:
                df = df.append(pd.DataFrame({'HZ_pt':df_[key_]['ak8']['H_pt'], 'NN':df_[key_]['val']['NN'], 'Weight':get_weight(cut), 'Name':name}), ignore_index=True)
            return(df)
        #
        if (sig in key_):
            ttZbb   = (df_[key_]['val']['Zbb'] == True)
            ttHbb   = (df_[key_]['val']['Hbb'] == True)
            ttZqq   = (df_[key_]['val']['Zqq'] == True)
            #
            sig_df = add_toDF(sig_df,ttZbb,'ttZbb')
            sig_df = add_toDF(sig_df,ttHbb,'ttHbb')
            bkg_df = add_toDF(bkg_df,ttZqq,'ttZqq')
        else:
            bkg_df = add_toDF(bkg_df)
    #
    # Seperate NN histogram by H/Z pt bin
    pt_bins = [200,300] # [200,300,400]
    def binSigBkg(df_):
        hist_vals  = []
        for i,pt_bin in enumerate(pt_bins):
            if i == len(pt_bins)-1:
                pt_cut = df_['HZ_pt'] > pt_bin
            else:
                pt_cut = (df_['HZ_pt'] >= pt_bins[i]) & (df_['HZ_pt'] < pt_bins[i+1])
            temp_vals,_ = np.histogram(df_['NN'][pt_cut], bins=50, range=(0,1), weights=df_['Weight'][pt_cut])
            hist_vals.append(temp_vals)
        return np.array(hist_vals)
    sigH_vals = binSigBkg(sig_df[sig_df['Name'] == 'ttHbb'])
    sigZ_vals = binSigBkg(sig_df[sig_df['Name'] == 'ttZbb'])
    bkg_vals  = binSigBkg(bkg_df)
    #
    # Calculate signal strength per bin
    def calcMu(sig_,bkg_):
        def GLS(data_,x_,bkg_=None):
            #
            X_     = np.matrix(x_)
            Omega_ = np.matrix(np.diag(data_))
            Y_     = np.matrix(data_-bkg_).T
            #
            Cov_   = (X_ * Omega_.I * X_.T).I # transpose is flipped because X is inheirently flipped in shape
            beta_  = Cov_ * (X_ * Omega_.I * Y_)
            Corr_  = Cov_.A / np.sqrt(np.outer(np.diag(Cov_),np.diag(Cov_)))
            return beta_.A.flatten(),Cov_,Corr_
        #
        mu_binned     = []
        for s_, b_, in zip(sig_,bkg_):
            mu_     = []
            for _ in range(10000):
                pData = np.random.poisson(s_+b_)
                x_,cov_,corr_ = GLS(pData,s_,b_)
                mu_.append(x_)
            mu_binned.append(mu_)
        return np.array(mu_binned).T
    #
    def plotMu(mu_,zorh):
        plt.hist(mu_, histtype='step', range=(-5,5), label=['200t300','gt300']) # ['200t300','300t400','gt400']
        plt.title(zorh)
        plt.legend()
        plt.show()
        plt.clf()
        #
        mu     = np.median(mu_, axis = 0)
        sigma1 = np.percentile(mu_,[16,84], axis=0)
        sigma2 = np.percentile(mu_,[5,95], axis=0)
        x_bins, x_bin_errs = [250,400], [50,100] #[250,350,500], [50,50,100]
        plt.errorbar(x=x_bins,y=mu,
                     xerr=x_bin_errs,
                     yerr=[-1*(sigma2[0]-mu),sigma2[1]-mu],
                     fmt='o',barsabove=True,capsize=5, label='2 sigma')
        plt.errorbar(x=x_bins,y=mu,
                     xerr=x_bin_errs,
                     yerr=[-1*(sigma1[0]-mu),sigma1[1]-mu],
                     fmt='o',barsabove=True,capsize=5, label='1 sigma')
        plt.title('tt'+zorh+', H/Z to bb')
        plt.xlabel(zorh+' pt (GeV)')
        plt.ylabel('$\mu$')
        plt.legend()
        # plt.ylim(0,None)
        plt.xlim(200,500)
        plt.show()
    #
    HZ_sigstr = np.squeeze(calcMu(sigH_vals+sigZ_vals, bkg_vals))
    Z_sigstr = np.squeeze(calcMu(sigZ_vals,bkg_vals+sigH_vals))
    H_sigstr = np.squeeze(calcMu(sigH_vals,bkg_vals+sigZ_vals))
    plotMu(HZ_sigstr,'Z/H')
    plotMu(Z_sigstr ,'Z')
    plotMu(H_sigstr, 'H')
    #
def diffXSec(files_, samples_, outDir_):
    df_ = kFit.retrieveData(files_, samples_, outDir_, getak8_ = True)
    sig = 'TTZH'
    suf = '_2017'
    #
    sig_df = pd.DataFrame()
    bkg_df = pd.DataFrame()
    for key_ in df_.keys():
        def get_weight(cut):
            weight = (df_[key_]['val']['weight']*np.sign(df_[key_]['val']['genWeight']) * (137/41.9))[cut]
            return weight
        def add_toDF(df, cut, name=key_.split('_201')[0]):
            df = df.append(pd.DataFrame({'HZ_pt':df_[key_]['ak8']['H_pt'][cut], 'Weight':get_weight(cut), 'Name':name}), ignore_index=True)
            return(df)
        #
        base_cuts =(
            (df_[key_]['ak8']['n_nonHbb'] >= 2)    &
            (df_[key_]['ak8']['nhbbFatJets'] > 0)  &
            (df_[key_]['val']['NN']         >= .96)&
            (df_[key_]['ak8']['H_pt']       >= 200)&
            (df_[key_]['ak8']['H_M']         > 50) &
            (df_[key_]['ak8']['H_M']         < 200))
        
        if (sig in key_):
            ttZbb   = base_cuts & (df_[key_]['val']['Zbb'] == True)
            ttHbb   = base_cuts & (df_[key_]['val']['Hbb'] == True)
            ttZqq   = base_cuts & (df_[key_]['val']['Zqq'] == True)
            #
            sig_df = add_toDF(sig_df,ttZbb,'ttZbb')
            sig_df = add_toDF(sig_df,ttHbb,'ttHbb')
            bkg_df = add_toDF(bkg_df,ttZqq,'ttZqq')
        else:
            bkg_df = add_toDF(bkg_df,base_cuts)
    #
    bins = np.array([200,300,500])
    def binSigBkg(arr_,w_, raw_=False):
        arr_ = np.clip(arr_,bins[0],bins[-1])
        hist_vals,_ = np.histogram(arr_,bins=bins,weights=w_)
        if raw_:
            raw_vals,_ = np.histogram(arr_,bins=bins)
            return hist_vals,raw_vals
        return hist_vals
    #
    sigH_vals, rawSigH = binSigBkg(sig_df['HZ_pt'][sig_df['Name'] == 'ttHbb'],sig_df['Weight'][sig_df['Name'] == 'ttHbb'],True)
    sigZ_vals, rawSigZ = binSigBkg(sig_df['HZ_pt'][sig_df['Name'] == 'ttZbb'],sig_df['Weight'][sig_df['Name'] == 'ttZbb'],True)
    bkg_vals = binSigBkg(bkg_df['HZ_pt'],bkg_df['Weight'])  
    print('\n',rawSigH)
    print(rawSigZ,'\n')
    #
    #
    def plotSignalStrength(xsec_,sig_eff,title_,extra_bkg=0):
        ttZH_200to300 = []
        ttZH_gt300    = []
        for _ in range(10000):
            pdata_vals = np.random.poisson(sigH_vals+sigZ_vals+bkg_vals)
            pSig = pdata_vals - bkg_vals - extra_bkg
            Xsec_ttzh = pSig/137/1000/sig_eff
            ttZH_200to300.append(Xsec_ttzh[0]/xsec_)
            ttZH_gt300.append(Xsec_ttzh[1]/xsec_)
        #
        plt.hist(ttZH_200to300,histtype='step',label='200 to 300')
        plt.hist(ttZH_gt300,histtype='step',label='>= 300')
        plt.title(title_)
        plt.xlabel('$\mu$')
        plt.legend()
        plt.show()
        plt.clf()
        #
        plt.errorbar(x=[250,400],y=np.array([np.median(ttZH_200to300),np.median(ttZH_gt300)]),
                     xerr=[50,100],yerr=np.array([2*np.std(ttZH_200to300),2*np.std(ttZH_gt300)]),
                     fmt='o',barsabove=True,capsize=5, label='2 sigma')
        plt.errorbar(x=[250,400],y=np.array([np.median(ttZH_200to300),np.median(ttZH_gt300)]),
                     xerr=[50,100],yerr=np.array([np.std(ttZH_200to300),np.std(ttZH_gt300)]),
                     fmt='o',barsabove=True,capsize=5, label='1 sigma')
        plt.title(title_)
        plt.xlabel('Z/H pt (GeV)')
        plt.ylabel('$\mu$')
        plt.legend()
        # plt.ylim(0,None)
        plt.xlim(200,500)
        plt.show()
    #
    plotSignalStrength(cfg.ZHbbtotXsec,(rawSigH + rawSigZ)/cfg.n_ZHbbMC,'ttZ/H, Z/H to bb')
    plotSignalStrength(cfg.ZHbbXsec['ttZbb'], rawSigZ/cfg.n_ZHbbMC_dict['ttZbb'], 'ttZ, Z to bb', extra_bkg = sigH_vals)
    plotSignalStrength(cfg.ZHbbXsec['ttHbb'], rawSigH/cfg.n_ZHbbMC_dict['ttHbb'], 'ttH, H to bb', extra_bkg = sigZ_vals)

if __name__ == '__main__':
    files_samples_outDir = cfg.ZHbbFitCfg
    #
    #Bkg_Est(*files_samples_outDir, cfg.ZHbbFitoverlap)
    measSigStr(*files_samples_outDir)
    #diffXSec(*files_samples_outDir)
