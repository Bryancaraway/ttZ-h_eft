import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import re
from lib.fun_library import save_pdf
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg

#pull_txt = f'{sys.path[1]}/test/fitdiag_roots/partblind_pull_4.txt'
#name = 'partblind_uncorrbtaglfhf_pull_'
name = 'partblind_pull_final'
#name = 'partblind_pull_alt'
#name = 'partblind_pull_ouZH_b1_pt'
#name = 'partblind_pull_Zh_M_NNcuts'
pull_txt = f'{sys.path[1]}/test/fitdiag_roots/{name}.txt'

def main():
    pulls = {}
    with open(pull_txt, 'r') as pull_file:
        # parse name, pull, +/-
        for l in pull_file.readlines():
            if "name" in l or '[1.00, 1.00]' in l: continue # skip first line
            #print(l.split()) # ['CMS_ttbbnorm', '[0.00,', '5.00]', '+1.37', '+/-', '0.20', '+1.37', '+/-', '0.20', '+0.00']
            if 'CMS' in l:
                name, nom, bands = l.split()[0], l.split()[3], l.split()[5]
            else:
                name, nom, bands = l.split()[0], l.split()[4], l.split()[6]
            #print(name,nom,bands)
            pulls[name] = [float(nom.strip('*!')),float(bands)]
            
    pull_df = pd.DataFrame.from_dict(pulls,orient='index',columns=['nom','bands'])
    #print(pull_df)
    #
    make_nuicance_pull_plots(pull_df)
    
@save_pdf(f'nuicance_{name}.pdf')
def make_nuicance_pull_plots(df):
    # first get jes and btag
    jec_btagw_df = pd.concat([
        df.loc[df.index.str.contains('jes')],
        df.loc[df.index.str.contains('jer')],
        df.loc[df.index.str.contains('jms')],
        df.loc[df.index.str.contains('jmr')],
        df.loc[df.index.str.contains('btg')],
    ])
    stat_df = df.loc[df.index.str.contains('prop')]
    rest_df = df.drop(jec_btagw_df.index).drop(stat_df.index)
    jec_btagw_df = jec_btagw_df[~jec_btagw_df.index.duplicated(keep='last')]
    #
    fig, ax = initplt()
    make_errorbar(fig, ax, jec_btagw_df)
    fig, ax = initplt()
    make_errorbar(fig, ax, rest_df)
    fig, ax = initplt()
    make_errorbar(fig, ax, stat_df)


def make_errorbar(fig,ax,df_):
    ax.errorbar(
        x=np.arange(0,len(df_),1), y=df_['nom'], 
        yerr=df_['bands'],
        fmt='.', color='k', label='Fit to Sig+Bkg (blinded)')
    ax.tick_params(which='both', direction='in')
    ax.set_xticks(np.arange(0,len(df_),1))
    ax.set_xticklabels(df_.index, rotation=90, fontsize=4)
    endplt(fig,ax, len(df_))

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
    ax.fill_between([-1,300], y1=1, y2=-1, color='k', alpha=0.1)
    ax.set_xlim(-1,max_)
    ax.set_ylim(-2.5,2.5)
    ax.set_ylabel(r'$\theta$')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis='x', color='k', linestyle=':', alpha=0.25)
    ax.legend()


if __name__ == '__main__':
    main()
