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

stxs_exp = '''
   r_ttZ0 :    +1.000   -9.467/+9.505 (68%)
   r_ttZ1 :    +1.000   -2.124/+2.260 (68%)
   r_ttZ2 :    +1.000   -1.326/+1.394 (68%)
   r_ttZ3 :    +1.000   -1.730/+1.929 (68%)
   r_ttH0 :    +1.000   -5.852/+5.653 (68%)
   r_ttH1 :    +1.000   -7.787/+8.265 (68%)
   r_ttH2 :    +1.000   -1.259/+1.329 (68%)
   r_ttH3 :    +1.000   -1.359/+1.505 (68%)
'''
# this hasnt been updated 
stxs_300_exp = '''
   r_ttZ0 :    +1.000   -7.476/+7.695 (68%)
   r_ttZ1 :    +1.000   -1.990/+2.097 (68%)
   r_ttZ2 :    +1.000   -0.924/+0.990 (68%)
   r_ttH0 :    +1.000   -5.571/+5.425 (68%)
   r_ttH1 :    +1.000   -7.097/+7.344 (68%)
   r_ttH2 :    +1.000   -0.797/+0.846 (68%)
'''

inc_exp = '''
   r_ttZ :    +1.000   -0.861/+0.936 (68%)
   r_ttH :    +1.000   -0.730/+0.789 (68%)
'''


def plot_exp_xsec(p_xs, p_exp_mu, bins, p_type):
    print(p_xs)

    y_stxs      = np.array([p_xs[p_type+str(i)]['yield'] for i in range(1,4)])
    y_stxs_up   = np.array([p_xs[p_type+str(i)]['yield']*p_exp_mu[p_type+str(i)][0] for i in range(1,4)])
    y_stxs_down = np.array([p_xs[p_type+str(i)]['yield']*p_exp_mu[p_type+str(i)][1] for i in range(1,4)])
    # get theory uncertainties from json
    y_theo_up   = np.array([p_xs[p_type+str(i)]['theo_Up']-p_xs[p_type+str(i)]['yield'] for i in range(1,4)])
    y_theo_down = np.array([p_xs[p_type+str(i)]['yield'] - p_xs[p_type+str(i)]['theo_Down'] for i in range(1,4)])
    
    #y_inc       = np.array([sum([p_xs[p_type+str(i)] for i in range(1,4) ])])
    #y_inc_up    = np.array([sum([p_xs[p_type+str(i)] for i in range(1,4) ]) * p_exp_mu[p_type][0]])
    #y_inc_down  = np.array([sum([p_xs[p_type+str(i)] for i in range(1,4) ]) * p_exp_mu[p_type][1]])
    fig, ax = beginPlt()
    #ax.errorbar(
    #    x=(bins['stxs'][1:]+bins['stxs'][:-1])/2, 
    #    y=y_stxs,
    #    #xerr=(bins['stxs'][1:]-bins['stxs'][:-1])/2, 
    #    yerr=[abs(y_stxs_down),y_stxs_up],
    #    fmt='.', color='r', label='STXS'
    #)
    ax.step(
        x=bins['stxs'],
        y=np.append(y_stxs,y_stxs[-1]),
        where='post',color='k',label='Prediction'
    )
    make_error_boxes(ax, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_stxs, 
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, np.array([abs(y_stxs_down),y_stxs_up]),
                     facecolor='gray', alpha=0.3, hatch=None, label='Pred unc.')
    #make_error_boxes(ax, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_stxs, 
    #                 (bins['stxs'][1:]-bins['stxs'][:-1])/2, np.array([abs(y_theo_down),y_theo_up]),
    #                 facecolor='green', alpha=0.3, hatch=None, label='Theo unc.')
    #ax.set_xticks(bins['stxs'])
    endPlt(ax,p_type, ylabel=r'$\sigma$ [fb]', errorbox=True)
    plt.tight_layout()
    #plt.show()
    #
    #fig, ax = beginPlt()
    #y_stxs = np.append(y_stxs,0)
    #y_inc = np.append(y_inc,0)
    #print(y_stxs)
    #print(y_inc)
    #y_stxs = abs(y_stxs[1:]-y_stxs[:-1])/(bins['stxs'][1:]-bins['stxs'][:-1])
    #y_inc = abs(y_inc[1:]-y_inc[:-1])/(bins['inc'][1:]-bins['inc'][:-1])
    #print(y_stxs)
    #print(y_inc)
    #exit()
    

def plot_exp_mu(mu_str, bin_dict, llabel):
    mu_dict = {}
    for l in mu_str.split('\n'):
        if len(l) == 0: continue
        mu_dict.update(parse_mu_line(l))
    #
    fig, ax = beginPlt()
    #for i,mu in enumerate(mu_dict):
    mu_up   = np.array([mu_dict[mu][0] for mu in mu_dict])
    mu_down = np.array([mu_dict[mu][1] for mu in mu_dict])
    ax.errorbar(x=np.ones(len(mu_dict)), 
                xerr=[abs(mu_down),mu_up], y=np.arange(len(mu_dict)), 
                fmt='.', c='r',label=llabel)
    ax.axvline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
    plt.legend()
    ax.set_yticks(np.arange(len(mu_dict)))
    ax.set_yticklabels(y_label_format(mu_dict, bin_dict), fontsize=10)
    ax.set_xticks([-2,-1,0,1,2,3,4])
    ax.set_xlabel(r'best-fit $\mu$')
    ax.set_xlim(-2.2, 4.2)
    ax.set_ylim(-1, len(mu_dict))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(axis='x')
    #plt.show()

def y_label_format(mu_dict, bin_dict):
    l_list = []
    for mu in mu_dict:
        bin_num = re.search(r'(\d)?$',mu).group()
        l1 = f"{bin_dict[bin_num]} GeV\n"
        p = re.search(r'tt(Z|H)',mu).group()
        l2 = rf"${{\mu}}_{{\mathrm{{{p}}}}} = 1.00_{{{mu_dict[mu][1]}}}^{{+{mu_dict[mu][0]}}}$"
        l_list.append(l1+l2)
    return l_list

def beginPlt():
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.88,wspace=0.0)
    #fig.text(0.15,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 10)
    #fig.text(0.70,0.89, f'137'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
    CMSlabel(fig,ax)
    return fig,ax

def endPlt(ax,p_type, ylabel, doLog=True, errorbox=False):
    #
    ax.set_ylabel(ylabel)
    process = re.search(r'(Z|H)',p_type).group()
    ax.set_xlabel(rf"Simulated ${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{process}}}}}$ [GeV]", usetex=True)
    #
    ax.tick_params(which='both', direction='in', top=True, right=True)
    #
    ax.set_xlim(200,600)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xticklabels(f'{t:.0f}' if t != 600 else r'$\infty$' for t in ax.get_xticks())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if doLog:
        ax.set_yscale('Log')
    ax.set_ylim(1 if ax.get_yscale() == 'log' else 0, ax.get_ylim()[1] * (10 if ax.get_yscale() == 'log' else 1.5))
    if errorbox:
        handles, labels = ax.get_legend_handles_labels()
        asi_patch = Patch(label='Asimov Data',  fc='gray', snap=True, alpha=0.3)
        theo_patch = Patch(label='Theo unc.',  fc='green', snap=True, alpha=0.3)
        #handles = [(asi_patch,handles[0]),(theo_patch,handles[0])]#[(asi_patch,handles[0])]
        handles = [(asi_patch,handles[0])]
        labels  = labels#+['Expected']
        ax.legend(handles,labels, #bbox_to_anchor=(1.00,1), 
                  fontsize=10, framealpha = 0, loc='upper right')
    else:
        plt.legend()
    #plt.show()

@save_pdf("exp_diffxsec.pdf")
def make_exp_xsec(p_xs, p_exp_mu):
    bins = {'inc':np.array([200,600]),'stxs':np.array([200,300,450,600])}
    for p in p_xs:
        plot_exp_xsec(p_xs[p], p_exp_mu, bins, p)

def get_p_xs():
    _out_dict = {}
    p_df = json.load(open(cfg.dataDir+'/process_norms/process_norms_ttbbw_run2.json','r'))
    # also get mu_rf, mu_r, mu_f, pdf, alphas
    scale = ['mu_f','mu_r','mu_rf']
    pdfas = ['pdfweight','alphas']
    for p in ['ttZ','ttH']:
        o_df = {**{p:0},**{p+str(i):{'yield':0, 'theo_Up':0, 'theo_Down':0} for i in range(4)}}
        for y in cfg.Years:
            o_df[p] = o_df[p]+p_df[y][p]['stxs_yield']
            for i in range(4):
                o_df[p+str(i)]['yield'] = o_df[p+str(i)]['yield']+p_df[y][p+str(i)]['stxs_yield']
                for ud_str in ['Up','Down']:
                    scale_unc = max([abs(1-1/p_df[y][p+str(i)][f'{sc}_{ud_str}']) for sc in scale])
                    pdfas_unc = np.sqrt((1-1/p_df[y][p+str(i)][f'pdfweight_{ud_str}'])**2 + (1-1/p_df[y][p+str(i)][f'alphas_{ud_str}'])**2)
                    theo_unc = np.sqrt((scale_unc)**2 + (pdfas_unc)**2)
                    if ud_str == 'Up':
                        theo_unc = 1+theo_unc
                    else:
                        theo_unc = 1-theo_unc
                    o_df[p+str(i)][f'theo_{ud_str}'] = o_df[p+str(i)][f'theo_{ud_str}']+p_df[y][p+str(i)]['stxs_yield']*theo_unc
            #
        #
        xs_df = {**{f'{p}':signal_xsec[p]*1000}, **{f'{p}{i}': {k:signal_xsec[p]*o_df[p+str(i)][k]/o_df[p]*1000 for k in o_df[p+str(i)]} for i in range(4)}}
        _out_dict[p] = xs_df
    return _out_dict

def parse_p_exp_mu():
    _out_dict = {}
    for l in stxs_exp.split('\n'):
        if len(l) == 0: continue
        _out_dict.update(parse_mu_line(l))
    for l in inc_exp.split('\n'):
        if len(l) == 0: continue
        _out_dict.update(parse_mu_line(l))
    return _out_dict
    
def parse_mu_line(l):
    p = re.search(r'tt(Z|H)\d*',l).group()
    ud  = l.split()[3].split('/')
    return {p:[float(ud[1]),float(ud[0])]}

@save_pdf('exp_mu.pdf')
def make_exp_mu():
    plot_exp_mu(inc_exp, {'':"Inc"}, 'Inc')
    plot_exp_mu(stxs_300_exp, {'0':"[0, 200]", '1':"[200, 300]", '2':"[300, inf]"}, 'STXS')
    plot_exp_mu(stxs_exp, {'0':"[0, 200]", '1':"[200, 300]", '2':"[300, inf]", '3':"[450, inf]"}, 'STXS')
           
def main():
    import_mpl_settings(1)
    p_xs = get_p_xs()
    p_exp_mu = parse_p_exp_mu()
    #
    make_exp_xsec(p_xs, p_exp_mu)
    make_exp_mu()

if __name__ == '__main__':
    main()

