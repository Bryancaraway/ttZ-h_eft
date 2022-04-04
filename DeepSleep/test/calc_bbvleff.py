####  CALC LEPTON bbvl eff/fake ####
### Written by: Bryan Caraway    ###
####################################

import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
#import functools
import config.ana_cff as cfg
from config.sample_cff import sample_cfg, process_cfg
from lib.fun_library import clop_pear_ci, save_pdf, import_mpl_settings, CMSlabel


sigs = ['TTZToBB','ttHTobb']
bkgs    = ['QCD_HT_500to700','QCD_HT_700to1000','QCD_HT_1000to1500','QCD_HT_1500to2000','QCD_HT_2000toInf']
#pt_bins = np.array([200, 300, 400, 500, 700, 1000])
pt_bins = np.array([200, 300, 450, 1000])

@save_pdf('bbvl_eff_mistag.pdf')
def handle_bbvl_eff_mistag():
    for year in cfg.Years:
        for sig in sigs:
            s_info = get_data(sig,year)
            eff = calc_bbvl_eff_mistag(s_info)
            print(sig, 'INC', sum(s_info['tot_num'])/sum(s_info['tot_den']))
            print(sig, '200-450', sum(s_info['tot_num'][0:2])/sum(s_info['tot_den'][0:2]))
            print(sig, '>450', s_info['tot_num'][-1]/s_info['tot_den'][-1],'\n')
            plot_eff_mistag(eff, sig, year, 'eff')
        for bkg in bkgs:
            b_info = get_data(bkg,year)
            mistag = calc_bbvl_eff_mistag(b_info)
            plot_eff_mistag(mistag, bkg, year, 'mistag')

def plot_eff_mistag(i_dict, sample, year, i_type):
    for cat in ['msd','bbvl','both']:
        if 'QCD' in sample and cat =='msd':
            continue
        fig, ax = beginPlt(year)
        ax.errorbar(x=(pt_bins[1:]+pt_bins[:-1])/2, y=i_dict[cat][0],
                     xerr=(pt_bins[1:]-pt_bins[:-1])/2, 
                     yerr=[i_dict[cat][0] - i_dict[cat][1], i_dict[cat][2] - i_dict[cat][0]],#down,up
                     fmt='.', color='k',label=f'{sample} {year}')
        ax.set_ylabel(f'{cat} {i_type}')
        ax.set_xlabel(r'AK8 jet $p_{\mathrm{T}}$ [GeV]')
        ax.set_ylim((0,1.1) if i_type == 'eff' else (0,.1))
        endPlt(ax)

def beginPlt(y):
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.88,wspace=0.0,hspace=0.0)
    CMSlabel(fig,ax,opt='Simulation', lumi=round(cfg.Lumi[y],1))
    return fig, ax
def endPlt(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xlim(0,1000)
    ax.legend()
    plt.tight_layout()

def calc_bbvl_eff_mistag(_df):
    _out = {
        'msd' :[_df['msd_num']/_df['tot_den'], *clop_pear_ci(_df['msd_num'] ,_df['tot_den'])],
        'bbvl':[_df['bbvl_num']/_df['tot_den'],*clop_pear_ci(_df['bbvl_num'],_df['tot_den'])],
        'both':[_df['tot_num']/_df['tot_den'], *clop_pear_ci(_df['tot_num'] ,_df['tot_den'])],
    }
    return _out

def get_data(_s,_y):
    return pd.read_pickle(cfg.postSkim_dir+f"{_y}/bbvL_{sample_cfg[_s]['out_name']}/{_s}.pkl")['metaData']

if __name__ == '__main__':
    import_mpl_settings()
    handle_bbvl_eff_mistag()
