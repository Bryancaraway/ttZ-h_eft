import uproot
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
#mass = 'Zh_M_alt'
mass = 'Zh_M'

#nn = 'NN' # just to test

target = [mass,'Zh_pt']
doNNcuts = True
tbins_map = {
    mass : np.arange(50,200+5,5), # new format
    #'Zh_M' :[50,75,90,105,120,140,200], # new format
    'Zh_pt':[200,300,450,np.inf],
    'Zh_score':[0,.25,.50,1.0],
    'nBottoms_drLeptonCleaned':[1.5,2.5,3.5,4.5,10],
    'n_ak4jets':[3.5,4.5,5.5,6.5,7.5,14],
    'n_b_inZh':[-0.5,0.5,1.5,2.5]}

def initDF(data, add_cut=(lambda _df: _df['Zh_pt'] >= 0), add_cut_str=''):
    #processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','VJets','other']
    #processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','single_t','ttX','VJets']
    processes = ['ttZ','ttH','TTBar','tt_B','single_t','ttX','VJets']
    processes += ['data_obs'] # adding data
    #rare_p = {'ttX','VJets','other'}    
    k_list = [mass,'Zh_pt',nn,'process','n_ak4jets']
    df = pd.DataFrame()
    for p in processes:
        for y in ['2016','2017','2018']:
            if p != 'data_obs':
                data[f'{p}_{y}']['event_weight'] = getZhbbWeight(data[f'{p}_{y}'],y)
                #if p == 'tt_B':
                #    data[f'{p}_{y}']['event_weight'] = 1.5*data[f'{p}_{y}']['event_weight']
                data[f'{p}_{y}']['event_weight2'] = data[f'{p}_{y}']['event_weight']**2

                    
            else: # is data
                data[f'{p}_{y}']['event_weight'] = np.ones(len(data[f'{p}_{y}']['Zh_pt']))
                data[f'{p}_{y}']['event_weight2'] = data[f'{p}_{y}']['event_weight']**2
            #if p in rare_p:
            #    data[f'{p}_{y}']['process'] = 'rare'
            df = pd.concat([df,data[f'{p}_{y}'].filter(items=k_list+['event_weight','event_weight2'])], 
                           axis='rows', ignore_index=True)
    #
    # organize by process
    #cuts = (lambda df_: df_[(df_[nn] > 0.80) & (df_['Zh_pt'] > pt_cut) & (df_['n_ak4jets'] >= 5)])
    cuts = (lambda df_: df_[(df_[nn] > 0.80) & (add_cut(df_)) & (df_['n_ak4jets'] >= 5)])
    df = cuts(df)
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
        l,c   = np.array([np.array(getLaLabel(p, altcolors=True)) for p in p_]).T
        if issig:
            c = [_.replace('tab:orange','magenta') for _ in c]
        #l   = np.array([f'{x} ({y:3.1f}+/-{z:3.1f})' for x,y,z in zip(l,i,i_err)])
        if issig is not None:
            #xfactor = [int(issig/i_p) for i_p in i]
            xfactor = [10 for i_p in i]
            w = np.array([w_p*x_p for w_p,x_p in zip(w,xfactor)], dtype=object)
            l   = np.array([rf'{l_p}$\times{x_p}$' for l_p,x_p in zip(l,xfactor)])
        return h,w,w2,i,c,l
    #
    
    bkg_h, bkg_w, bkg_w2, bkg_i, bkg_c, bkg_l = get_hist_essentials(bkg_p)
    
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
    n, _ = np.histogram(np.hstack(bkg_h),bins=bins, weights=np.hstack(bkg_w),range=(bins[0],bins[-1]))
    x,xerr,y,yerr =getMCStat_info(n, np.hstack(bkg_h),np.hstack(bkg_w), bins)
    make_error_boxes(ax,x,y,xerr,yerr,label='Stat unc.')
    #plot step sig
    sig_h, sig_w, sig_w2, sig_i, sig_c, sig_l = get_hist_essentials(sig_p,issig=sum(n))
    sig_ls = [':','--']
    for i_ in range(len(sig_h)):
        _ = ax.hist(
            sig_h[i_],
            bins=tbins_map[mass],
            histtype='step',
            linewidth = 1.0, # i think default is 1
            linestyle = sig_ls[i_],
            weights   = sig_w[i_],
            color     = sig_c[i_],
            label     = sig_l[i_]
    )
    #
    #plot step data
    # ============
    #data_n, edges = np.histogram(df[mass][df['process'] == 'Data'].values, bins = bins, range=(bins[0],bins[-1]))
    #ax.errorbar(
    #    x=(edges[1:] + edges[:-1])/2 , y=data_n,
    #    xerr=(edges[1:] - edges[:-1])/2, yerr=np.sqrt(data_n),
    #    fmt='.', label='data', color='k')
    #
    endplt(fig,ax,add_cut_str)
    #plt.show()
    #plt.savefig(f'pdf/pas_zhm_80nn_{pt_cut}pt_run2.pdf')


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

def endplt(fig,ax,add_cut_str):
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    #
    CMSlabel(fig=fig, ax=ax, opt='Simulation', lumi='nl')
    #CMSlabel(fig=fig, ax=ax, opt='Simulation Preliminary', lumi='nl')
    tex_x_corr = 0.55 if '> \mathsf{450}' not in add_cut_str else 0.62
    tex_x_corr = 0.53 if '\mathsf{200} <' in add_cut_str else tex_x_corr
    #fig.text(tex_x_corr,0.59, rf'{{{add_cut_str}}} GeV', usetex=True, fontsize=6)
    #fig.text(tex_x_corr,0.55, r'DNN score $>\mathsf{0.80}$', usetex=True, fontsize=6)
    fig.text(tex_x_corr,0.59, rf'{{{add_cut_str}}} GeV', usetex=True, fontsize=7)
    fig.text(tex_x_corr,0.53, r'DNN score $>\mathsf{0.80}$', usetex=True, fontsize=7)
    #ax.set_xlabel(r'${m}_{\mathrm{SD}}^{\mathrm{Z/H\; cand.}}$ [GeV]', usetex=True)
    ax.set_xlabel(r'$\mathsf{m}_{\text{SD}}^{\text{Z/H\;cand.}}$ \raisebox{0.25ex}{[}$\text{GeV}$\raisebox{0.25ex}{]}', usetex=True)
    #self.ax.set_ylabel(f"{'%' if self.doNorm else 'Events'} / {(self.bin_w[0].round(2) if len(np.unique(self.bin_w.round(4))) == 1 else 'bin')}")#fontsize = self.fontsize)
    ax.set_ylabel('Events / 5 GeV', usetex=True) # hardcoded
    #print(self.ax.get_ylim()[1], self.ax.get_ylim()[1] * 1.10 )        
    #plt.xlim(self.bin_range)
    #ax.set_yscale('log')
    ax.set_xlim(tbins_map[mass][0],tbins_map[mass][-1])
    #
    y_scale_upper = 1.70 if '\mathsf{300} <' in add_cut_str else 1.5
    y_scale_upper = 1.2 if '\mathsf{200} <' in add_cut_str else y_scale_upper
    ax.set_ylim(0,ymax=ax.get_ylim()[1]*y_scale_upper)
    #plt.grid(True)
    handles, labels = ax.get_legend_handles_labels()
    hatch_patch = Patch(hatch=10*'X', label='Stat Unc.',  fc='w', alpha=0.99)
    handles = handles + [hatch_patch]
    labels  = labels + ['Stat Unc.']
    ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize=7)
    plt.tight_layout()
    #plt.show()

#@save_pdf(f'pas_zhm_80newnn_ptcut_run2.pdf')
@save_pdf(f'pas_zhm_80newnn_ptcut_run2_final.pdf')
def main():
    pas_data_file = cfg.dataDir+'/pas_plot_info/pas_data_file.pkl'
    if not os.path.exists(pas_data_file):
        data = PlotsFromDC(
            sig= cfg.Sig_MC,
            bkg = cfg.Bkg_MC,
            #extrak= ['newgenm_NN'] + ([] if mass == 'Zh_M' else [mass]),
            extrak= [] if mass == 'Zh_M' else [mass],
            extrasigk=['matchedGen_ZHbb_bb']).getData()
        AnaDict(data).to_pickle(pas_data_file)
    else:
        data = AnaDict.read_pickle(pas_data_file)
    #
    import_mpl_settings(1)
    
    initDF(data, 
           add_cut=(lambda _df: ((_df['Zh_pt'] >= 200) & (_df['Zh_pt'] < 300)) ),
           #add_cut_str=r'$200 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 300$')
           add_cut_str=r'$\mathsf{200} < \mathsf{p}_{\text{T}}^{\text{Z/H\;cand.}} < \mathsf{300}$')
    initDF(data,
           add_cut=(lambda _df: ((_df['Zh_pt'] >= 300) & (_df['Zh_pt'] < 450)) ), 
           #add_cut_str=r'$300 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 450$')
           add_cut_str=r'$\mathsf{300} < \mathsf{p}_{\text{T}}^{\text{Z/H\;cand.}} < \mathsf{450}$')
    initDF(data,
           add_cut=(lambda _df: ((_df['Zh_pt'] >= 450) & (_df['Zh_pt'] < np.inf)) ), 
           #add_cut_str=r'${p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} > 450$')
           add_cut_str=r'$\mathsf{p}_{\text{T}}^{\text{Z/H\;cand.}} > \mathsf{450}$')


if __name__ == '__main__':
    main()
