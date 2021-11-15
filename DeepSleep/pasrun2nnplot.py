import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)',
        shell=True).decode().strip('\n')+'DeepSleep/')
#
import re
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

nn = cfg.nn
#nn = 'newgenm_NN'

def plot_nn(data, add_cut=(lambda _df: _df['Zh_pt'] >= 0), add_cut_str=''):
    processes = ['ttZ','ttH','TTBar','tt_B']
    sig_p = ['ttZ','ttH']
    bkg_p = ['TTBar','tt_B']
    k_list = ['Zh_pt',nn,'process','matchedGen_ZHbb_bb']
    df = pd.DataFrame()
    for p in processes:
        for y in ['2016','2017','2018']:
            data[f'{p}_{y}']['event_weight'] = getZhbbWeight(data[f'{p}_{y}'],y)
            data[f'{p}_{y}']['event_weight2'] = data[f'{p}_{y}']['event_weight']**2
            #if p in rare_p:
            #    data[f'{p}_{y}']['process'] = 'rare'
            df = pd.concat([df,data[f'{p}_{y}'].filter(items=k_list+['event_weight','event_weight2'])], 
                           axis='rows', ignore_index=True)

    cuts = (lambda df_: df_[add_cut(df_)])
    df = cuts(df)
    
    def get_hist_essentials(p_,issig=None):
        if issig:
            getv = (lambda k: 
                    [df[k][(df['process'] == p) & (df['matchedGen_ZHbb_bb'] == True)].values for p in p_] + \
                    [df[k][(df['process'] == p) & (df['matchedGen_ZHbb_bb'] == False)].values for p in p_] )
            l_p = [f'{p}_genm_'+re.search(r'(Z|H)',p).group()+'bb_bb' for p in p_] + [f'{p}_notgenm_'+re.search(r'(Z|H)',p).group()+'bb' for p in p_]
        else:
            getv = (lambda k: [df[k][df['process'] == p].values for p in p_])
            l_p = p_
        #
        h = getv(nn)
        w = getv('event_weight')
        w2 = getv('event_weight2')
        # done with getv
        i     = [np.sum(w_i) for w_i in w]
        i_err     = [np.sqrt(np.sum(w2_i)) for w2_i in w2]
        #w = np.divide(w,i)   # to nomralize
        #w2 = np.divide(w2,i) # to normalize
        #i_err = np.sqrt(np.sum(w2,axis=0))
        l,c   = np.array([np.array(getLaLabel(p)) for p in l_p]).T
        c = [_.replace('gold','magenta') for _ in c]
        #l   = np.array([f'{x} ({y:3.1f}+/-{z:3.1f})' for x,y,z in zip(l,i,i_err)])
        return h,w,w2,i,c,l
    #
    fig,ax = initplt()
    bins = np.linspace(0,1,21)
    def plotter(pro_h,pro_w,pro_w2,pro_i,pro_c,pro_l):
        for i_ in range(len(pro_h)):
            n_ = ax.hist(
                pro_h[i_],
                bins=bins,
                histtype='step',
                #linewidth = 1.5, # i think default is 1
                #linestyle = pro_ls[i_],
                weights   = pro_w[i_]/pro_i[i_],
                color     = pro_c[i_],
                label     = pro_l[i_]
            )
            y    = np.histogram(pro_h[i_], bins=bins, weights=pro_w[i_])[0]
            yerr = np.histogram(pro_h[i_], bins=bins, weights=pro_w2[i_])[0]
            ax.errorbar(
                x=(bins[1:]+bins[:-1])/2, y=y/pro_i[i_],
                yerr=np.sqrt(yerr)/pro_i[i_],
                fmt='.',ms=3, color=pro_c[i_]
            )
    #bkg_h, bkg_w, bkg_w2, bkg_i, bkg_c, bkg_l = get_hist_essentials(bkg_p)
    #sig_h, sig_w, sig_w2, sig_i, sig_c, sig_l = get_hist_essentials(sig_p, issig=True)
    plotter(*get_hist_essentials(sig_p, issig=True))
    plotter(*get_hist_essentials(bkg_p))
    #for i_ in range(len(bkg_h)):
    #    n_ = ax.hist(
    #        bkg_h[i_],
    #        bins=bins,
    #        histtype='step',
    #        #linewidth = 1.5, # i think default is 1
    #        #linestyle = bkg_ls[i_],
    #        weights   = bkg_w[i_],
    #        color     = bkg_c[i_],
    #        label     = bkg_l[i_]
    #    )
    #    y    = np.histogram(bkg_h[i_], bins=bins, weights=bkg_w[i_])[0]
    #    yerr = np.histogram(bkg_h[i_], bins=bins, weights=bkg_w2[i_])[0]
    #    ax.errorbar(
    #        x=(bins[1:]+bins[:-1])/2, y=y,
    #        yerr=np.sqrt(yerr),
    #        fmt='.',ms=3, color=bkg_c[i_]
    #    )

    endplt(fig,ax,add_cut_str)


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
    #self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Simulation$", fontsize = self.fontsize)
    #fig.text(0.105,0.89, r"\textbf{CMS} {\footnotesize \textit{Simulation}}", usetex=True, fontsize = 10)
    #fig.text(0.635,0.89, f'137'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
    CMSlabel(fig=fig, ax=ax, opt='Simulation')
    fig.text(0.50,0.63, rf'{{{add_cut_str}}} GeV', usetex=True, fontsize=6)
    #fig.text(0.635,0.62, r'DNN score $>0.80$', usetex=True, fontsize=10)
    ax.set_xlabel(r'DNN score', usetex=True)
    #self.ax.set_ylabel(f"{'%' if self.doNorm else 'Events'} / {(self.bin_w[0].round(2) if len(np.unique(self.bin_w.round(4))) == 1 else 'bin')}")#fontsize = self.fontsize)
    ax.set_ylabel('fraction of yield / bin', usetex=True) # hardcoded
    #print(self.ax.get_ylim()[1], self.ax.get_ylim()[1] * 1.10 )        
    #plt.xlim(self.bin_range)
    ax.set_yscale('log')
    #ax.set_xlim(tbins_map['Zh_M'][0],tbins_map['Zh_M'][-1])
    ax.set_xlim(0,1)
    #ax.set_ylim(ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1]*7.5)
    ax.set_ylim(0.001,ymax=ax.get_ylim()[1]*15)
    #plt.grid(True)
    ax.legend(framealpha = 0, ncol=2, fontsize=6)
    plt.tight_layout()
    #plt.show()

@save_pdf(f'pas_nncomp_ptcut_run2.pdf')   
def main():
    pas_data_file = cfg.dataDir+'/pas_plot_info/pas_data_file.pkl'
    data = AnaDict.read_pickle(pas_data_file)
    import_mpl_settings(1) 

    plot_nn(
        data,
        add_cut=(lambda _df: ((_df['Zh_pt'] >= 400) & (_df['Zh_pt'] < np.inf)) ),
        add_cut_str=r'${p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} > 300$'
    )

if __name__ == '__main__':
    main()
