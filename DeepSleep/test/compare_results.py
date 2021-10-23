import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

import json
import numpy as np
import pandas as pd
import re

import config.ana_cff as cfg
from config.sample_cff import signal_xsec
from lib.fun_library import t2Run, save_pdf, getZhbbBaseCuts, getZhbbWeight, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel, make_error_boxes


pt_bins = [0,200,300,450,600]


inc_exp_2016 = '''
   r_ttZ :    +1.000   -1.000/+1.628 (68%)
   r_ttH :    +1.000   -1.000/+1.513 (68%)
'''
inc_exp_2017 = '''
    r_ttZ :    +1.000   -1.000/+1.949 (68%)
    r_ttH :    +1.000   -1.000/+1.594 (68%)
'''
inc_exp_2018 = '''
    r_ttZ :    +1.000   -1.000/+1.377 (68%)
    r_ttH :    +1.000   -1.000/+1.130 (68%)
'''
inc_exp_run2 = '''
    r_ttZ :    +1.000   -0.866/+0.939 (68%)
    r_ttH :    +1.000   -0.734/+0.796 (68%)
'''
inc_obs_2016 = '''
   r_ttZ :    +2.703   -1.775/+2.093 (68%)
   r_ttH :    +0.463   -0.463/+1.746 (68%)
'''
inc_obs_2017 = '''
    r_ttZ :    +0.786   -0.786/+2.359 (68%)
    r_ttH :    +0.423   -0.423/+1.844 (68%)
'''
inc_obs_2018 = '''
    r_ttZ :    +0.000   -0.000/+0.679 (68%)
    r_ttH :    +0.000   -0.000/+0.549 (68%)
'''
inc_obs_run2 = '''
    r_ttZ :    +0.502   -0.502/+1.058 (68%)
    r_ttH :    +0.000   -0.000/+0.643 (68%)
'''
stxs_exp_2016 = '''
   r_ttZ0 :    +1.000   -1.000/+20.439 (68%)
   r_ttZ1 :    +1.000   -1.000/+4.019 (68%)
   r_ttZ2 :    +1.000   -1.000/+2.416 (68%)
   r_ttZ3 :    +1.000   -1.000/+4.717 (68%)
   r_ttH0 :    +1.000   -1.000/+7.983 (68%)
   r_ttH1 :    +1.000   -1.000/+14.001 (68%)
   r_ttH2 :    +1.000   -1.000/+2.278 (68%)
   r_ttH3 :    +1.000   -1.000/+3.031 (68%)
'''
stxs_exp_2017 = '''
   r_ttZ0 :    +1.000   -1.000/+14.499 (68%)
   r_ttZ1 :    +1.000   -1.000/+4.537 (68%)
   r_ttZ2 :    +1.000   -1.000/+3.201 (68%)
   r_ttZ3 :    +1.000   -1.000/+3.520 (68%)
   r_ttH0 :    +1.000   -1.000/+8.275 (68%)
   r_ttH1 :    +1.000   -1.000/+13.571 (68%)
   r_ttH2 :    +1.000   -1.000/+2.384 (68%)
   r_ttH3 :    +1.000   -1.000/+3.020 (68%)
'''
stxs_exp_2018 = '''
   r_ttZ0 :    +1.000   -1.000/+13.200 (68%)
   r_ttZ1 :    +1.000   -1.000/+2.983 (68%)
   r_ttZ2 :    +1.000   -1.000/+2.055 (68%)
   r_ttZ3 :    +1.000   -1.000/+2.830 (68%)
   r_ttH0 :    +1.000   -1.000/+7.288 (68%)
   r_ttH1 :    +1.000   -1.000/+11.348 (68%)
   r_ttH2 :    +1.000   -1.000/+1.695 (68%)
   r_ttH3 :    +1.000   -1.000/+2.253 (68%)
'''
stxs_exp_run2 = '''
   r_ttZ0 :    +1.000   -1.000/+9.285 (68%)
   r_ttZ1 :    +1.000   -1.000/+2.217 (68%)
   r_ttZ2 :    +1.000   -1.000/+1.394 (68%)
   r_ttZ3 :    +1.000   -1.000/+1.929 (68%)
   r_ttH0 :    +1.000   -1.000/+5.055 (68%)
   r_ttH1 :    +1.000   -1.000/+7.731 (68%)
   r_ttH2 :    +1.000   -1.000/+1.196 (68%)
   r_ttH3 :    +1.000   -1.000/+1.505 (68%)
'''
stxs_obs_2016 = '''
   r_ttZ0 :    +9.416   -9.416/+20.584 (68%)
   r_ttZ1 :    +1.971   -1.971/+4.819 (68%)
   r_ttZ2 :    +4.770   -2.892/+3.348 (68%)
   r_ttZ3 :    +0.000   -0.000/+2.510 (68%)
   r_ttH0 :    +0.419   -0.419/+10.814 (68%)
   r_ttH1 :    +1.613   -1.613/+17.997 (68%)
   r_ttH2 :    +2.349   -2.349/+2.707 (68%)
   r_ttH3 :    +0.000   -0.000/+1.066 (68%)
'''
stxs_obs_2017 = '''
   r_ttZ0 :    +0.002   -0.002/+13.283 (68%)
   r_ttZ1 :    +0.000   -0.000/+1.694 (68%)
   r_ttZ2 :    +3.111   -3.111/+4.362 (68%)
   r_ttZ3 :    +2.010   -2.010/+4.176 (68%)
   r_ttH0 :    +0.000   -0.000/+3.020 (68%)
   r_ttH1 :    +0.000   -0.000/+2.825 (68%)
   r_ttH2 :    +2.072   -2.072/+2.778 (68%)
   r_ttH3 :    +0.000   -0.000/+2.202 (68%)
'''
stxs_obs_2018 = '''
   r_ttZ0 :    +0.000   -0.000/+9.239 (68%)
   r_ttZ1 :    +0.000   -0.000/+1.570 (68%)
   r_ttZ2 :    +0.000   -0.000/+1.877 (68%)
   r_ttZ3 :    +0.000   -0.000/+1.414 (68%)
   r_ttH0 :    +0.012   -0.012/+6.538 (68%)
   r_ttH1 :    +0.000   -0.000/+4.921 (68%)
   r_ttH2 :    +0.000   -0.000/+1.365 (68%)
   r_ttH3 :    +0.000   -0.000/+0.842 (68%)
'''
stxs_obs_run2 = '''
   r_ttZ0 :    +0.000   -0.000/+9.019 (68%)
   r_ttZ1 :    +0.000   -0.000/+1.031 (68%)
   r_ttZ2 :    +2.029   -1.541/+1.628 (68%)
   r_ttZ3 :    +0.000   -0.000/+1.492 (68%)
   r_ttH0 :    +0.000   -0.000/+2.143 (68%)
   r_ttH1 :    +0.000   -0.000/+1.933 (68%)
   r_ttH2 :    +0.978   -0.978/+1.281 (68%)
   r_ttH3 :    +0.000   -0.000/+0.455 (68%)=
'''
# unconstrained
#inc_exp_2016 = '''
#    r_ttZ :    +1.000   -1.451/+1.629 (68%)
#    r_ttH :    +1.000   -1.333/+1.513 (68%)
#'''
#inc_exp_2017 = '''
#    r_ttZ :    +1.000   -1.791/+1.949 (68%)
#    r_ttH :    +1.000   -1.423/+1.594 (68%)
#'''
#inc_exp_2018 = '''
#    r_ttZ :    +1.000   -1.273/+1.377 (68%)
#    r_ttH :    +1.000   -1.045/+1.130 (68%)
#'''
#inc_exp_run2 = '''
#    r_ttZ :    +1.000   -0.866/+0.939 (68%)
#    r_ttH :    +1.000   -0.734/+0.796 (68%)
#'''
#inc_obs_2016 = '''
#    r_ttZ :    +2.628   -1.881/+2.118 (68%)
#    r_ttH :    +0.386   -1.646/+1.787 (68%)
#'''
#inc_obs_2017 = '''
#    r_ttZ :    -0.339   -2.798/+2.684 (68%)
#    r_ttH :    -0.407   -2.045/+1.983 (68%)
#'''
#inc_obs_2018 = '''
#    r_ttZ :    -1.157   -1.459/+1.433 (68%)
#    r_ttH :    -0.951   -1.246/+1.248 (68%)
#'''
#inc_obs_run2 = '''
#    r_ttZ :    +0.214   -1.052/+1.091 (68%)
#    r_ttH :    -0.492   -0.890/+0.900 (68%)
#'''
#
#stxs_exp_2016 = '''
#
#'''
#stxs_exp_2017 = '''
#'''
#stxs_exp_2018 = '''
#'''
stxs_exp_run2 = '''
   r_ttZ0 :    +1.000   -9.466/+9.503 (68%)
   r_ttZ1 :    +1.000   -2.124/+2.260 (68%)
   r_ttZ2 :    +1.000   -1.326/+1.394 (68%)
   r_ttZ3 :    +1.000   -1.730/+1.929 (68%)
   r_ttH0 :    +1.000   -5.852/+5.653 (68%)
   r_ttH1 :    +1.000   -7.786/+8.265 (68%)
   r_ttH2 :    +1.000   -1.259/+1.329 (68%)
   r_ttH3 :    +1.000   -1.359/+1.505 (68%)
'''
#stxs_obs_2016 = '''
#   r_ttZ0 :    +7.106   -33.391/+22.894 (68%)
#   r_ttZ1 :    +1.920   -4.601/+5.082 (68%)
#   r_ttZ2 :    +4.926   -3.053/+3.379 (68%)
#   r_ttZ3 :    -1.189   -3.456/+4.159 (68%)
#   r_ttH0 :    -0.438   -11.680/+11.111 (68%)
#   r_ttH1 :    +0.053   -16.828/+18.625 (68%)
#   r_ttH2 :    +2.809   -2.908/+3.200 (68%)
#   r_ttH3 :    -3.247   -1.919/+2.507 (68%)
#'''
#stxs_obs_2017 = '''
#
#'''
#stxs_obs_2018 = '''
#   r_ttZ0 :   -17.602   -12.398/+23.423 (68%)
#   r_ttZ1 :    -2.842   -3.409/+3.498 (68%)
#   r_ttZ2 :    +0.069   -2.259/+2.294 (68%)
#   r_ttZ3 :    -1.060   -2.646/+2.869 (68%)
#   r_ttH0 :   +11.538   -13.693/+14.595 (68%)
#   r_ttH1 :   -13.907   -16.093/+17.370 (68%)
#   r_ttH2 :    +0.952   -2.289/+2.369 (68%)
#   r_ttH3 :    -2.377   -2.062/+2.329 (68%)
#'''
stxs_obs_run2 = '''
   r_ttZ0 :   +12.466   -16.216/+15.259 (68%)
   r_ttZ1 :    +0.231   -2.598/+2.726 (68%)
   r_ttZ2 :    +1.208   -1.737/+1.799 (68%)
   r_ttZ3 :    -0.237   -1.894/+2.050 (68%)
   r_ttH0 :    -2.488   -8.875/+8.955 (68%)
   r_ttH1 :   -20.793   -9.207/+10.908 (68%)
   r_ttH2 :    +3.042   -1.619/+1.719 (68%)
   r_ttH3 :    -3.163   -1.382/+1.473 (68%)
'''

@save_pdf("compare_results_run2_unconstrained.pdf")
def main():
    bin_dict = {'0':"[0, 200]", '1':"[200, 300]", '2':"[300, 450]", '3':r"[450, \inf]"}
    inc_bin_dict = {'':"Inc"}
    #plot_exp_mu(stxs_exp, bin_dict,       'murf uncorrelated (current)', fig, ax, 0)
    #plot_exp_mu(xtrasmooth_exp, bin_dict, 'Proposed (rate uncertainty)', fig, ax, 1)
    #plot_exp_mu(murf_exp, bin_dict, 'murf correlated', fig, ax, 1)
    #
    #fig, ax = beginPlt('2016')
    #plot_exp_mu(inc_exp_2016, inc_bin_dict,       'Expected', fig, ax, 0)  
    #plot_exp_mu(inc_obs_2016, inc_bin_dict,       'Observed', fig, ax, 1)  
    #plt.legend(fontsize=8)
    ##
    #fig, ax = beginPlt('2016')
    #plot_exp_mu(stxs_exp_2016, bin_dict,       'Expected', fig, ax, 0)
    #plot_exp_mu(stxs_obs_2016, bin_dict,       'Observed', fig, ax, 1)
    #plt.legend(fontsize=8)
    ##
    #fig, ax = beginPlt('2017')
    #plot_exp_mu(inc_exp_2017, inc_bin_dict,       'Expected', fig, ax, 0)  
    #plot_exp_mu(inc_obs_2017, inc_bin_dict,       'Observed', fig, ax, 1)  
    #plt.legend(fontsize=8)
    ##
    #fig, ax = beginPlt('2017')
    #plot_exp_mu(stxs_exp_2017, bin_dict,       'Expected', fig, ax, 0)
    #plot_exp_mu(stxs_obs_2017, bin_dict,       'Observed', fig, ax, 1)
    #plt.legend(fontsize=8)
    ##
    #fig, ax = beginPlt('2018')
    #plot_exp_mu(inc_exp_2018, inc_bin_dict,       'Expected', fig, ax, 0)  
    #plot_exp_mu(inc_obs_2018, inc_bin_dict,       'Observed', fig, ax, 1)  
    #plt.legend(fontsize=8)
    ##
    #fig, ax = beginPlt('2018')
    #plot_exp_mu(stxs_exp_2018, bin_dict,       'Expected', fig, ax, 0)
    #plot_exp_mu(stxs_obs_2018, bin_dict,       'Observed', fig, ax, 1)
    #plt.legend(fontsize=8)
    ##
    #fig, ax = beginPlt()
    #plot_exp_mu(inc_exp_run2, inc_bin_dict,       'Expected', fig, ax, 0)  
    #plot_exp_mu(inc_obs_run2, inc_bin_dict,       'Observed', fig, ax, 1)  
    #plt.legend(fontsize=8)
    ##
    fig, ax = beginPlt()
    plot_exp_mu(stxs_exp_run2, bin_dict,       'Expected', fig, ax, 0)
    plot_exp_mu(stxs_obs_run2, bin_dict,       'Observed', fig, ax, 1)
    plt.legend(fontsize=8)
    #
    #plt.show()

def plot_exp_mu(mu_str, bin_dict, llabel, fig, ax, i):
    c_list = ['r','b']
    mu_dict = {}
    for l in mu_str.split('\n'):
        if len(l) == 0: continue
        mu_dict.update(parse_mu_line(l))
    #
    #for i,mu in enumerate(mu_dict):
    mu_nom  = np.array([mu_dict[mu][0] for mu in mu_dict])
    mu_up   = np.array([mu_dict[mu][1] for mu in mu_dict])
    mu_down = np.array([mu_dict[mu][2] for mu in mu_dict])
    ax.errorbar(x=mu_nom,
                #x=np.ones(len(mu_dict)), 
                xerr=[abs(mu_down),mu_up], y=np.arange(len(mu_dict))-(i*.1*np.ones(len(mu_dict))), 
                fmt='.', c=c_list[i],label=llabel)


    if i == 0:
        ax.axvline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        ax.axhline(3.5, color='k', linewidth='1', linestyle='-', snap=True)
        ax.set_yticks(np.arange(len(mu_dict)))
        ax.set_yticklabels(y_label_format(mu_dict, bin_dict), fontsize=10)
        ax.set_xlabel(r'best-fit $\mu$')
        #ax.set_xticks([-2,-1,0,1,2,3,4])
        #ax.set_xlim(-2, 4)
        ax.set_xticks(np.arange(-8,8,2))
        ax.set_xlim(-8,8)
        ax.set_ylim(-1, len(mu_dict))
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(axis='x')
    #plt.show()

def y_label_format(mu_dict, bin_dict):
    l_list = []
    for mu in mu_dict:
        bin_num = re.search(r'(\d)?$',mu).group()
        l1 = f"{bin_dict[bin_num]}"
        l1 += " GeV\n" if "Inc" not in l1 else " \n"
        p = re.search(r'tt(Z|H)',mu).group()
        l2 = rf"${{\mu}}_{{\mathrm{{{p}}}}}$"# = 1.00_{{{mu_dict[mu][1]}}}^{{+{mu_dict[mu][0]}}}$"
        l_list.append(l1+l2)
    return l_list

def parse_mu_line(l):
    p = re.search(r'tt(Z|H)\d*',l).group()
    nom = l.split()[2]
    ud  = l.split()[3].split('/')
    return {p:[float(nom),float(ud[1]),float(ud[0])]}

def beginPlt(year=None):
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.18,right=0.95,wspace=0.0)
    #fig.text(0.15,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 10)
    #fig.text(0.70,0.89, f'137'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
    CMSlabel(fig,ax,lumi=round(cfg.Lumi.get(year,137),1))
    return fig,ax

if __name__ == '__main__':
    import_mpl_settings(2)
    main()
