import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import re
from lib.fun_library import save_pdf, import_mpl_settings, CMSlabel
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg

#pull_txt = f'{sys.path[1]}/test/fitdiag_roots/partblind_pull_4.txt'
#name = 'partblind_pull_final'
name = 'unblind_pull_final'

pull_txt = f'{sys.path[1]}/test/fitdiag_roots/{name}.txt'

def main():
    pulls = {}
    with open(pull_txt, 'r') as pull_file:
        # parse name, pull, +/-
        for l in pull_file.readlines():
            if "name" in l or '[1.00, 1.00]' in l or 'r_tt' in l: continue # skip first line
            #print(l.split()) # ['CMS_ttbbnorm', '[0.00,', '5.00]', '+1.37', '+/-', '0.20', '+1.37', '+/-', '0.20', '+0.00']
            l_arr = l.split()
            name = l_arr[0]
            if 'CMS' in l:
                nom, bands = l_arr[6], float(l_arr[8])
            elif re.search(r"prop_biny\d*_Zhpt\d_bin\d*_\w*",l_arr[0]) and float(l_arr[1]) > 0: # poisson prior
                nom = str(1-(float(l_arr[9].strip('*!'))/float(l_arr[1]))) # to not break code
                up_band = float(l_arr[10])/float(l_arr[2])
                dn_band = abs(float(l_arr[11])/float(l_arr[3]))
                bands=[dn_band,up_band]
            else:
                nom, bands = l_arr[9], float(l_arr[11])
            nom = float(nom.strip('*!'))
            pulls[name] = [nom,bands]
            
    pull_df = pd.DataFrame.from_dict(pulls,orient='index',columns=['nom','bands'])
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
    #special handling for single value bands
    stat_poisson_df = df.filter(regex=r'\w*bin\d*_\w*',axis='index')
    stat_poisson_df.loc[:,'bands'] = stat_poisson_df['bands'].apply( 
        func=(lambda _x : list([_x,_x]) if type(_x) is float else _x)).copy()
    #
    stat_df = stat_df.drop(stat_poisson_df.index) # drop poisson uncertainties
    jec_btagw_df = jec_btagw_df[~jec_btagw_df.index.duplicated(keep='last')]
    #
    print(stat_df[abs(stat_df['nom']) > 0.7].index)
    #
    fig, ax = initplt()
    make_errorbar(fig, ax, jec_btagw_df)
    fig, ax = initplt()
    make_errorbar(fig, ax, rest_df)
    fig, ax = initplt()
    make_errorbar(fig, ax, stat_df)
    fig, ax = initplt()
    make_errorbar(fig, ax, stat_poisson_df)


def make_errorbar(fig,ax,df_):
    # handling for list yerr
    if np.array(df_['bands'].to_list()).ndim == 1:
        yerr = df_['bands']
    else: # handle asym errorbars
        yerr = np.array(df_['bands'].to_list())
        ynom = df_['nom'].to_numpy()
        yerr = (ynom+yerr.T)
    ax.errorbar(
        x=np.arange(0,len(df_),1), y=df_['nom'], 
        yerr=yerr,
        fmt='.', color='k', label='Fit to Sig+Bkg',zorder=2)#label='Fit to Sig+Bkg (blinded)')
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
    CMSlabel(fig, ax, fontsize=12)
    ax.axhline(0,  color='red', linewidth='1', linestyle='-', snap=True, alpha=0.5)
    ax.fill_between([-1,300], y1=1, y2=-1, color='k', alpha=0.1,zorder=1)
    ax.set_xlim(-1,max_)
    ax.set_ylim(-2.5,2.5)
    ax.set_ylabel(r'$\theta_\mathrm{NP}$')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis='x', color='k', linestyle=':', alpha=0.25)
    ax.legend()
    plt.tight_layout()


if __name__ == '__main__':
    import_mpl_settings(no_figsize=True)
    main()
