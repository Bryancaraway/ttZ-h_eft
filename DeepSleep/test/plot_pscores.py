import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
import numpy as np
np.random.seed(1)
import re
import uproot
from scipy.stats import norm
import config.ana_cff as cfg
from lib.fun_library import save_pdf, getLaLabel
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            

gof_dir = 'fitdiag_roots/gof'
fit_vars = cfg.withbbvl_dnn_ZHgenm_vars
fit_vars = [f'NN_{i}' for i in range(64)]

@save_pdf("qcnn_pscores_hl2.pdf")
def main():
    #dnn_vars = cfg.withbbvl_dnn_ZH_vars

    p_df = pd.DataFrame.from_dict({v:calc_p(v) for v in fit_vars},orient='index',columns=['pscore'])
    fig, ax = initplt()
    ax.errorbar(
        x=np.arange(0,len(p_df),1), y=p_df['pscore'], 
        #yerr=p_df['bands'],
        fmt='.', color='k', label='Fit to Sig+Bkg (blinded)')
    ax.tick_params(which='both', direction='in')
    ax.set_xticks(np.arange(0,len(p_df),1))
    ax.set_xticklabels(p_df.index, rotation=90, fontsize=4)
    endplt(fig,ax, len(p_df))
    plt.close()
    plt.hist(p_df['pscore'], bins=20, range=(0,1))
    for i in p_df['pscore']:
        print(i)
    #plt.show()
 
def calc_p(var):
    t     = uproot.open(f"{gof_dir}/higgsCombine_{var}_NNcuts_run2.GoodnessOfFit.mH120.root")['limit']
    t_toy = uproot.open(f"{gof_dir}/higgsCombine_{var}_NNcuts_run2_toy.GoodnessOfFit.mH120.root")['limit']
    gof_data =  t.array('limit')[0]
    gof_toy  =  t_toy.array('limit')
    pscore    = len(gof_toy[gof_toy>=gof_data])/len(gof_toy)
    return pscore

def initplt():
    fig, ax = plt.subplots()
    fig.subplots_adjust(
            top=0.88,
            bottom=0.20,
            left=0.11,
            right=0.88,
            hspace=0.2,
            wspace=0.2
    )
    return fig, ax

def endplt(fig,ax, max_):
    fig.text(0.105,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 12)
    fig.text(0.70,0.89, f'137'+r' fb$^{-1}$ (13 TeV)',  fontsize = 12)
    ax.axhline(0,  color='red', linewidth='1', linestyle='-', snap=True, alpha=0.5)
    ax.axhline(1,  color='red', linewidth='1', linestyle='-', snap=True, alpha=0.5)
    #ax.fill_between([-1,300], y1=1, y2=-1, color='k', alpha=0.1)
    ax.set_xlim(-1,max_)
    ax.set_ylim(-.20,1.20)
    ax.set_ylabel(r'p-score')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis='x', color='k', linestyle=':', alpha=0.25)
    ax.legend()



if __name__ == '__main__':
    main()
