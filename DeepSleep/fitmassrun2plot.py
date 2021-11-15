import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)',
        shell=True).decode().strip('\n')+'DeepSleep/')
#
import re
from multiprocessing import Pool
#from pathos.multiprocessing import ProcessingPool as Pool
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
import matplotlib.backends.backend_pdf as matpdf
from matplotlib import rc
rc("savefig",dpi=250)
rc("figure", figsize=(3.375, 3.375*(6./8.)), dpi=200)                                                            
rc("figure", max_open_warning=600)
rc("hatch", linewidth=0.5, color='r') 
#
import functools
import numpy as np
import pandas as pd
import config.ana_cff as cfg
from lib.fun_library import t2Run, save_pdf, getZhbbWeight, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel
from modules.AnaDict import AnaDict
from anrun2plot import PlotsFromDC

nn = cfg.nn
# testing random things ###
# fjetsdmnos_1 
#opt = 'nos'
#exopt = 'uncorr'
#opt  = 'alt'
##
#opt = ''
#exopt='corr'


# AN plots
opt='alt'
exopt=''
mass = f'fjetsdm{opt}_1'

#nn = 'NN' # just to test

target = [mass,'Zh_pt']
doNNcuts = True
tbins_map = {
    mass : np.linspace(60,120,61), # new format
 }

def initDF(data, add_cut=(lambda _df: _df['Zh_pt'] >= 0), add_cut_str='', years=cfg.Years):
    #processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','VJets','other']
    #processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','single_t','ttX','VJets']
    processes = ['ttZ','ttH','TTBar','tt_B','single_t','ttX','VJets','data_obs']
    #rare_p = {'ttX','VJets','other'}    
    k_list = [mass,'Zh_pt',nn,'process','n_ak4jets','fjetwmdscore_1']
    df = pd.DataFrame()
    for p in processes:
        for y in years:
            data[f'{p}_{y}']['event_weight'] = getZhbbWeight(data[f'{p}_{y}'],y) if p != 'data_obs' else np.ones_like(data[f'{p}_{y}'][mass])
            data[f'{p}_{y}']['event_weight2'] = data[f'{p}_{y}']['event_weight']**2 if p != 'data_obs' else np.ones_like(data[f'{p}_{y}'][mass])
            #if p in rare_p:
            #    data[f'{p}_{y}']['process'] = 'rare'
            df = pd.concat([df,data[f'{p}_{y}'].filter(items=k_list+['event_weight','event_weight2'])], 
                           axis='rows', ignore_index=True)
    #
    # organize by process
    #cuts = (lambda df_: df_[(df_[nn] > 0.80) & (df_['Zh_pt'] > pt_cut) & (df_['n_ak4jets'] >= 5)])
    cuts = (lambda df_: df_[add_cut(df_)])
    df = cuts(df)
    data_p = ['Data']
    sig_p = ['ttZ','ttH']
    #bkg_p = ['VJets','ttX','single_t','tt_2b','tt_bb','TTBar']
    bkg_p = ['VJets','ttX','single_t','tt_B','TTBar']
    #
    def get_hist_essentials(p_,issig=None):
        getv = (lambda k: [df[k][df['process'] == p].values for p in p_])
        #
        h = getv(mass)
        w = getv('event_weight')
        w2 = getv('event_weight2')
        i     = [np.sum(w_i) for w_i in w]
        i_err     = [np.sqrt(np.sum(w2_i)) for w2_i in w2]
        #i_err = np.sqrt(np.sum(w2,axis=0))
        l,c   = np.array([np.array(getLaLabel(p)) for p in p_]).T
        c = [_.replace('gold','magenta') for _ in c]
        #l   = np.array([f'{x} ({y:3.1f}+/-{z:3.1f})' for x,y,z in zip(l,i,i_err)])
        if issig is not None:
            xfactor = [int(issig/i_p) for i_p in i]
            w = np.array([w_p*x_p for w_p,x_p in zip(w,xfactor)], dtype=object)
            l   = np.array([rf'{l_p}$\times{x_p}$' for l_p,x_p in zip(l,xfactor)])
        return h,w,w2,i,c,l
    #
    from scipy.optimize import curve_fit
    def fit_func(x, l, q, x0, a, b, c): # fit function quadratic+gaussaian
        quadratic = q*np.power(x,2) + l*x + x0
        gaussian  = a*np.exp(-1*np.power(x-b,2)/(2*np.power(c,2)))
        return quadratic + gaussian

    #def fit_func(x, a, b, c): # fit function quadratic+gaussaian
    #    #quadratic = q*np.power(x,2) + l*x + x0
    #    gaussian  = a*np.exp(-1*np.power(x-b,2)/(2*np.power(c,2)))
    #    return gaussian
    

    bkg_h, bkg_w, bkg_w2, bkg_i, bkg_c, bkg_l = get_hist_essentials(sig_p+bkg_p)
    
    #
    fig,ax = initplt()
    #plot bkg stacked hist
    n_, bins, _ = ax.hist(
        bkg_h,
        bins=tbins_map[mass],
        stacked=True,
        histtype='stepfilled',
        weights = bkg_w,
        color   = bkg_c,
        label    = bkg_l
    )
    n_MC, _ = np.histogram(np.hstack(bkg_h),bins=bins, weights=np.hstack(bkg_w),range=(bins[0],bins[-1]))
    x,xerr,y,yerr =getMCStat_info(n_MC, np.hstack(bkg_h),np.hstack(bkg_w), bins)
    make_error_boxes(ax,x,y,xerr,yerr,label='Stat unc.')
    p0 = [(y[-1]-y[0])/(x[-1]-x[0]), -.09, 1, max(y) - 85*(y[-1]-y[0])/(x[-1]-x[0]), 85, 8]
    popt, pcov = curve_fit(fit_func, x, y, sigma=yerr, p0 = p0)
    popt_MC = popt
    x_sample = np.linspace(60,120,200)
    x_gaus = np.linspace(70,100,200)
    #print(popt)
    #ax.plot(x_sample, fit_func(x_sample, *p0), 'r-')
    #ax.plot(x_guas, 1200*np.exp(-1*np.power(x_gaus-85,2)/(2*np.power(15,2))), 'g-')
    ax.plot(x_sample, popt[1]*np.power(x_sample,2) + popt[0]*x_sample + popt[2], 'r-')
    ax.plot(x_sample, fit_func(x_sample, *popt), 'r-')
    fig.text(0.205,0.82, f'mean:{popt[4]:>10.2f}$\pm${pcov[4,4]**.5:.2f}\n$\sigma$:{abs(popt[5]):>10.3f}$\pm${pcov[5,5]**.5:.3f}', 
             usetex=True, fontsize=8, color='r')
    #
    data_h, data_w, data_w2, data_i, data_c, data_l = get_hist_essentials(data_p)
    n_data, _ = np.histogram(np.hstack(data_h),bins=bins, weights=np.hstack(data_w),range=(bins[0],bins[-1]))
    x,xerr,y,yerr = getMCStat_info(n_data, data_h, data_w, bins)
    ax.errorbar(x=x, y = y, xerr=xerr, yerr=yerr,
                fmt='.',  color=data_c[0], label=data_l[0])
    popt, pcov = curve_fit(fit_func, x, y, sigma=yerr, p0 = popt)
    popt_Data = popt
    #ax.plot(x_sample, fit_func(x_sample, *p0), 'b-')
    ax.plot(x_sample, popt[1]*np.power(x_sample,2) + popt[0]*x_sample + popt[2], 'b-')
    ax.plot(x_sample, fit_func(x_sample, *popt), 'b-')
    fig.text(0.205,0.74, f'mean:{popt[4]:>10.2f}$\pm${pcov[4,4]**.5:.2f}\n$\sigma$:{abs(popt[5]):>10.3f}$\pm${pcov[5,5]**.5:.3f}', 
             usetex=True, fontsize=8, color='b')
    #print(p0)
    #print(popt)
    print(popt_Data[4:6]/popt_MC[4:6])
    endplt(fig,ax,add_cut_str,years)

def getMCStat_info(n, h, w, bins):
    x = (bins[1:]+bins[:-1])/2
    xerr = (bins[1:] - bins[:-1])/2
    yerr = np.histogram(h, bins=bins, weights=np.power(w,2))[0]
    yerr = np.sqrt(yerr)
    y = n
    return x,xerr,y,yerr

def make_error_boxes(ax, xdata, ydata, xerror, yerror,  facecolor='r',
                         edgecolor='None', alpha=0.0, hatch=10*'X', label=''):
    errorboxes = []
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe, y - ye), 2*xe, 2*ye)
        errorboxes.append(rect)
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor, hatch=hatch, label=label, zorder=1.5)
    # Add collection to axes
    ax.add_collection(pc)


def initplt():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.2,
            wspace=0.2
        )
    return fig,ax

def endplt(fig,ax,add_cut_str,years=None):
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    CMSlabel(fig=fig, ax=ax, opt='Preliminary', lumi=round(cfg.Lumi[years[0]],1) if len(years) == 1 else None)
    #fig.text(0.635,0.66, rf'{{{add_cut_str}}} GeV', usetex=True, fontsize=10)
    fig.text(0.685,0.60, f'{add_cut_str} GeV', usetex=True, fontsize=8)
    #fig.text(0.635,0.62, r'DNN score $>0.80$', usetex=True, fontsize=10)
    ax.set_xlabel(rf'${{m}}_{{\mathrm{{SD({{{exopt}}})}}}}^{{\mathrm{{AK8\;jet}}}}$ [GeV]', fontsize = 10, usetex=True)
    ax.set_ylabel('Events / 5 GeV', fontsize=10, usetex=True) # hardcoded
    #plt.xlim(self.bin_range)
    ax.set_xlim(tbins_map[mass][0],tbins_map[mass][-1])
    #ax.set_ylim(ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1]*7.5)
    ax.set_ylim(0,ymax=ax.get_ylim()[1]*1.3)
    #plt.grid(True)
    handles, labels = ax.get_legend_handles_labels()
    hatch_patch = Patch(hatch=10*'X', label='Stat Unc.',  fc='w')
    handles = handles + [hatch_patch]
    labels  = labels + ['Stat Unc.']
    ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize=8, loc='upper right')
    plt.tight_layout()
    #plt.show()

@save_pdf(f'fitfjmass{opt}_run2.pdf')
def main():
    pas_data_file = cfg.dataDir+f'/pas_plot_info/fjmassfit{opt}_data_file.pkl'
    if not os.path.exists(pas_data_file):
        data = PlotsFromDC(
            sig= cfg.Sig_MC,
            bkg = cfg.Bkg_MC,
            cut = (lambda x_ : ((x_[mass] >= 60) & (x_[mass] <= 120) & 
                                (x_['fjetwmdscore_1'] > 0.0) & 
                                (x_['fjetpt_1'] > 200))),
            extrak= [mass,'fjetwmdscore_1','fjetpt_1'],
            extrasigk=['matchedGen_ZHbb_bb']).getData()
        AnaDict(data).to_pickle(pas_data_file)
    else:
        data = AnaDict.read_pickle(pas_data_file)
    #
    import_mpl_settings(1.5)
    #cut_str = r'$\mathrm{deepTagMD\;WvsQCD}>0.0$'+'\n'+r'${p}_{\mathrm{T}}^{\mathrm{AK8\;jet}} > 200$'
    cut_str = r'${p}_{\mathrm{T}}^{\mathrm{AK8\;jet}} > 200$'
    for year in cfg.Years:
        initDF(
            data, 
            #add_cut=(lambda x_ : (x_['fjetpt_1'] > 300)),
            add_cut = (lambda x_ : (x_['fjetwmdscore_1'] > 0.0)),
            add_cut_str=cut_str,
            years=[year],
        )
    initDF(
        data, 
        #add_cut=(lambda x_ : (x_['fjetpt_1'] > 300)),
        add_cut = (lambda x_ : (x_['fjetwmdscore_1'] > 0.0)),
        add_cut_str=cut_str,
        years=cfg.Years,
    )


if __name__ == '__main__':
    main()
