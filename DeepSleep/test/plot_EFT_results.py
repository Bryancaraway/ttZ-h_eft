import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc, lines
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

import numpy as np
import re

import config.ana_cff as cfg
from lib.fun_library import save_pdf, import_mpl_settings, upperlefttext, CMSlabel


wc_latex = {
    'cbW'  : r'$\frac{{c}_{\mathrm{bW}}}{{\Lambda}^{2}}$',
    'cptb' : r'$\frac{{c}_{\phi \mathrm{tb}}}{{\Lambda}^{2}}$',
    'cpt'  : r'$\frac{{c}_{\phi \mathrm{t}}}{{\Lambda}^{2}}$',
    'ctp'  : r'$\frac{{c}_{\mathrm{t} \phi}}{{\Lambda}^{2}}$',
    'ctZ'  : r'$\frac{{c}_{\mathrm{tZ}}}{{\Lambda}^{2}}$',
    'ctW'  : r'$\frac{{c}_{\mathrm{tW}}}{{\Lambda}^{2}}$',
    'cpQ3' : r'$\frac{{c}_{\phi \mathrm{Q}}^{3}}{{\Lambda}^{2}}$',
    'cpQM' : r'$\frac{{c}_{\phi \mathrm{Q}}^{-}}{{\Lambda}^{2}}$',
}
wc_bf_ffl = { # best-fit value
    'fixed':{
        'cbW' : -2.34,
        'cptb': 0.56,
        'cpt' : -0.26,
        'ctp' : 15.30,
        'ctZ' : 0.01,
        'ctW' : -0.05,
        'cpQ3': -0.54,
        'cpQM': -0.04,
    },
    'float':{
        'cbW' : 2.46,
        'cptb': 1.08,
        'cpt' : -0.34,
        'ctp' : 15.13,
        'ctZ' : -0.06,
        'ctW' : -0.07,
        'cpQ3': -0.27,
        'cpQM': 0.02,
    },
}
wc_seight_ffl = { # 68% CL lower, upper with respect to best-fit
    'fixed':{
        'cbW' : [-1.21,5.86],
        'cptb': [-5.46,5.28],
        'cpt' : [-3.54,2.99],
        'ctp' : [-8.65,8.50],
        'ctZ' : [-0.58,0.58],
        'ctW' : [-0.53,0.54],
        'cpQ3': [-1.77,1.80],
        'cpQM': [-2.47,2.70],
    },
    'float':{
        'cbW' : [-6.05,1.19],
        'cptb': [-5.63,5.25],
        'cpt' : [-5.22,4.39],
        'ctp' : [-8.15,8.08],
        'ctZ' : [-0.88,0.92],
        'ctW' : [-0.85,0.87],
        'cpQ3': [-2.00,2.03],
        'cpQM': [-3.66,4.02],
    },
}
wc_nfive_ffl = { # 95% CL lower, upper with respect to best-fit
    'fixed':{
        'cbW' : [-2.21,6.86],
        'cptb': [-10.00,9.59],
        'cpt' : [-8.24,5.65],
        'ctp' : [-15.17,14.88],
        'ctZ' : [-1.00,1.01],
        'ctW' : [-0.98,0.99],
        'cpQ3': [-3.32,3.42],
        'cpQM': [-4.73,5.67],
    },
    'float':{
        'cbW' : [-7.06,2.16],
        'cptb': [-10.53,9.58],
        'cpt' : [-10.73,7.91],
        'ctp' : [-14.63,14.50],
        'ctZ' : [-1.47,1.54],
        'ctW' : [-1.49,1.53],
        'cpQ3': [-3.77,3.83],
        'cpQM': [-6.71,7.95],
    },
}
@save_pdf('eft_results_table.pdf')
def main():
    fig, ax = beginPlt()
    for i,k in enumerate(wc_bf_ffl):
        plot_eft_result(ax,k,i)
    #
    endPlt(fig,ax)
    #plt.show()
    

def plot_eft_result(ax, f_fl, i_off):
    wc_seight = np.array([np.array(wc_seight_ffl[f_fl][wc]) + wc_bf_ffl[f_fl][wc] for wc in wc_latex])
    wc_nfive  = np.array([np.array(wc_nfive_ffl[f_fl][wc])  + wc_bf_ffl[f_fl][wc] for wc in wc_latex])
    y = np.arange(len(wc_latex))-(i_off*.3*np.ones(len(wc_latex)))+.15
    # plot 95% CL lines
    ax.hlines(
        y=y, xmin=wc_nfive[:,0], xmax=wc_nfive[:,1], 
        color=('r' if f_fl == 'fixed' else 'k'), linestyle='-', linewidth=.5)
    # plot 68% CL lines
    ax.hlines(
        y=y, xmin=wc_seight[:,0], xmax=wc_seight[:,1], 
        color=('r' if f_fl == 'fixed' else 'k'), linestyle='-', linewidth=2)
    

def beginPlt(year=None):
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.18,right=0.95,wspace=0.0)
    CMSlabel(fig,ax,lumi=round(cfg.Lumi.get(year,137),1))
    return fig,ax

def endPlt(fig,ax):
    ax.tick_params(which='both', direction='in', top=True)
    # yaxis
    ax.set_yticks(np.arange(len(wc_latex)))
    ax.set_yticklabels([wc_latex[wc] for wc in wc_latex], fontsize=16, usetex=True)
    # xaxis
    ax.set_xlabel(r'Wilson coefficient limit [$\mathrm{TeV}^{-2}$]', usetex=True)
    ax.set_xticks([-10, -5, 0, 5, 10, 15, 20, 25, 30])
    ax.set_xticklabels(['-10', '', '0', '', '10', '', '20', '', '30'])
    ax.set_xlim(-15,35)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis='x', linestyle=':')
    # legend
    handles = [
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='k'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='r'),
        lines.Line2D([],[], linestyle='-', linewidth=2, color='k'),
        lines.Line2D([],[], linestyle='-', linewidth=2, color='r'),
    ]
    labels  = [
        'Others profiled (95% CL)', 
        'Others fixed to SM (95% CL)', 
        'Others profiled (68% CL)', 
        'Others fixed to SM (68% CL)'
    ]
    ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.0, 0.72), fontsize=7, framealpha=1)

        

if __name__ == '__main__':
    import_mpl_settings(1, length=2, disable_sansmath=True)
    main()
