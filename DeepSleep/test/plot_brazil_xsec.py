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
brazil_green  = '#78c72c'
brazil_yellow = '#f9c754'

asym_fro_limits = { # r_ttZ0 and r_ttH0 float are fixed to SM
    'r_ttZ1':'''
    Observed Limit: r_ttZ1 < 3.3800
    Expected  2.5%: r_ttZ1 < 2.4596
    Expected 16.0%: r_ttZ1 < 3.3035
    Expected 50.0%: r_ttZ1 < 4.6641
    Expected 84.0%: r_ttZ1 < 6.6905
    Expected 97.5%: r_ttZ1 < 9.2952''',
    'r_ttZ2':'''
    Observed Limit: r_ttZ2 < 4.8860
    Expected  2.5%: r_ttZ2 < 1.6727
    Expected 16.0%: r_ttZ2 < 2.2466
    Expected 50.0%: r_ttZ2 < 3.1719
    Expected 84.0%: r_ttZ2 < 4.5247
    Expected 97.5%: r_ttZ2 < 6.2449''',
    'r_ttZ3':'''
    Observed Limit: r_ttZ3 < 4.0174
    Expected  2.5%: r_ttZ3 < 2.1350
    Expected 16.0%: r_ttZ3 < 2.8871
    Expected 50.0%: r_ttZ3 < 4.1406
    Expected 84.0%: r_ttZ3 < 6.0387
    Expected 97.5%: r_ttZ3 < 8.4702''',
    'r_ttH1':'''
    Observed Limit: r_ttH1 < 8.0311
    Expected  2.5%: r_ttH1 < 7.4956
    Expected 16.0%: r_ttH1 < 10.0274
    Expected 50.0%: r_ttH1 < 14.1094
    Expected 84.0%: r_ttH1 < 19.7897
    Expected 97.5%: r_ttH1 < 26.7422''',
    'r_ttH2':'''
    Observed Limit: r_ttH2 < 3.2512
    Expected  2.5%: r_ttH2 < 1.3290
    Expected 16.0%: r_ttH2 < 1.7923
    Expected 50.0%: r_ttH2 < 2.5391
    Expected 84.0%: r_ttH2 < 3.6422
    Expected 97.5%: r_ttH2 < 5.0283''',
    'r_ttH3':'''
    Observed Limit: r_ttH3 < 1.9635
    Expected  2.5%: r_ttH3 < 1.7000
    Expected 16.0%: r_ttH3 < 2.3113
    Expected 50.0%: r_ttH3 < 3.2969
    Expected 84.0%: r_ttH3 < 4.7556
    Expected 97.5%: r_ttH3 < 6.6288''',
}

asym_limits = asym_fro_limits

def plot_exp_xsec(p_xs, p_asymlim, bins, p_type):
    pt_ranges = np.array([[200,300],[300,450],[450,600]])
    #
    y_stxs      = np.array([p_xs[p_type+str(i)]['yield'] for i in range(1,4)])
    y_obs_lim   = np.array([p_xs[p_type+str(i)]['yield']*p_asymlim[p_type+str(i)]['obs'] for i in range(1,4)])
    y_exp_lim = np.array([p_xs[p_type+str(i)]['yield']*p_asymlim[p_type+str(i)]['exp'] for i in range(1,4)]) # 0 : -2sigma, 1: -1sigma, 2: nom, 3: +1sigma, 4: +2sigma
    print_latex_format(y_obs_lim, y_exp_lim, y_stxs, p_type)
    #print(p_type)
    #print("Median Expected")
    #print(y_exp_lim[:,2])
    #print(y_exp_lim[:,2]/y_stxs)
    #print("Median Expected + 1sigma")
    #print(y_exp_lim[:,-2]-y_exp_lim[:,2])
    #print((y_exp_lim[:,-2]-y_exp_lim[:,2])/y_stxs)
    #print("Median Expected - 1sigma")
    #print(y_exp_lim[:,2]-y_exp_lim[:,1])
    #print((y_exp_lim[:,2]-y_exp_lim[:,1])/y_stxs)
    # get theory uncertainties from json
    y_theo_up   = np.array([p_xs[p_type+str(i)]['theo_Up']-p_xs[p_type+str(i)]['yield'] for i in range(1,4)])
    y_theo_down = np.array([p_xs[p_type+str(i)]['yield'] - p_xs[p_type+str(i)]['theo_Down'] for i in range(1,4)])
    #
    fig, axs = beginPlt()
    ax = axs if len(axs) == '1' else axs[0]
    # Expected
    ax.hlines(y=y_exp_lim[:,2], xmin=pt_ranges[:,0], xmax=pt_ranges[:,1], color='#046cf4', linestyle='--', linewidth=0.5, label='Expected')
    make_error_boxes(ax, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_exp_lim[:,2], 
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, 
                     np.c_[y_exp_lim[:,2]-y_exp_lim[:,0],y_exp_lim[:,-1]-y_exp_lim[:,2]].T,
                     facecolor=brazil_yellow, alpha=1, hatch=None, label='Expected 95% CL')
    make_error_boxes(ax, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_exp_lim[:,2], 
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, 
                     np.c_[y_exp_lim[:,2]-y_exp_lim[:,1],y_exp_lim[:,-2]-y_exp_lim[:,2]].T,
                     facecolor=brazil_green, alpha=1, hatch=None, label='Expected 68% CL')

    # Observation
    ax.hlines(y=y_obs_lim, xmin=pt_ranges[:,0], xmax=pt_ranges[:,1], color='k', linestyle='-', linewidth=1, label='Observed')
    # Prediction
    ax.hlines(y=y_stxs, xmin=pt_ranges[:,0], xmax=pt_ranges[:,1], color='#ea7dfd',linestyle='-', linewidth=.7, label='SM')
    make_error_boxes(ax, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_stxs, 
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, np.array([abs(y_theo_down),y_theo_up]),
                     facecolor='#ea7dfd', alpha=0.3, hatch=None, label='SM cross-section')
    # Limit / Prediction
    ax_r = axs[1]
    # Observed
    ax_r.hlines(y=y_obs_lim/y_stxs, 
              xmin=pt_ranges[:,0], xmax=pt_ranges[:,1], color='k', linestyle='-', linewidth=1, label='Observed') 
    # Prediction
    ax_r.axhline(y=1, linestyle='-', linewidth=.7, color='#ea7dfd')
    make_error_boxes(ax_r, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_stxs/y_stxs, 
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, np.array([abs(y_theo_down/y_stxs),y_theo_up/y_stxs]),
                     facecolor='#ea7dfd', alpha=0.3, hatch=None, label='SM cross-section')
    # Expected
    ax_r.hlines(y=y_exp_lim[:,2]/y_stxs, xmin=pt_ranges[:,0], xmax=pt_ranges[:,1], color='#046cf4', linestyle='--', linewidth=0.5, label='Expected')
    make_error_boxes(ax_r, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_exp_lim[:,2]/y_stxs, 
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, 
                     (np.c_[y_exp_lim[:,2]-y_exp_lim[:,0],y_exp_lim[:,-1]-y_exp_lim[:,2]].T)/y_stxs,
                     facecolor=brazil_yellow, alpha=1, hatch=None, label='Expected 95% CL')
    make_error_boxes(ax_r, (bins['stxs'][1:]+bins['stxs'][:-1])/2, y_exp_lim[:,2]/y_stxs,  
                     (bins['stxs'][1:]-bins['stxs'][:-1])/2, 
                     (np.c_[y_exp_lim[:,2]-y_exp_lim[:,1],y_exp_lim[:,-2]-y_exp_lim[:,2]].T)/y_stxs,
                     facecolor=brazil_green, alpha=1, hatch=None, label='Expected 68% CL')
    #
    ax_r.set_ylabel('Limit / SM')
    #ax_r.yaxis.set_minor_locator(AutoMinorLocator())
    ax_r.tick_params(which='both', direction='in', top=True, right=True)
    ax_r.set_yscale('Log')
    ax_r.set_ylim(.5,19)
    #
    ax.set_xticks(bins['stxs'])
    #endPlt(ax,p_type, ylabel=r'$\sigma$ [fb]', errorbox=True) # no limit to SM ratio
    #endPlt(ax,p_type, ylabel=r'$\sigma$ [fb]', errorbox=True, ax_r=ax_r)
    endPlt(ax,p_type, ylabel=r'$\mathsf{\sigma}$ \raisebox{0.25ex}{[}$\text{fb}$\raisebox{0.25ex}{]}', errorbox=True, ax_r=ax_r)
    plt.tight_layout(h_pad=0)
    #plt.show()
    #exit()
    

def y_label_format(mu_dict, bin_dict):
    l_list = []
    for mu in mu_dict:
        bin_num = re.search(r'(\d)?$',mu).group()
        l1 = f"{bin_dict[bin_num]} GeV\n"
        p = re.search(r'tt(Z|H)',mu).group()
        #l2 = rf"${{\mu}}_{{\mathrm{{{p}}}}} = 1.00_{{{mu_dict[mu][1]}}}^{{+{mu_dict[mu][0]}}}$"
        l2 = rf"${{\mu}}_{{\mathsf{{{p}}}}} = 1.00_{{{mu_dict[mu][1]}}}^{{+{mu_dict[mu][0]}}}$"
        l_list.append(l1+l2)
    return l_list

def beginPlt():
    #fig, ax = plt.subplots()
    fig, axs = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[2,1]})
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.88,wspace=0.0,hspace=0.0)
    ax = axs if len(axs) == '1' else axs[0]
    #
    CMSlabel(fig,ax, opt='') # final, no preliminary
    #CMSlabel(fig,ax) 
    return fig, axs

def endPlt(ax,p_type, ylabel, doLog=True, errorbox=False, ax_r=None):
    #
    ax.set_ylabel(ylabel)
    process = re.search(r'(Z|H)',p_type).group()

    #
    ax.tick_params(which='both', direction='in', top=True, right=True)
    #
    ax.set_xlim(200,600)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if ax_r is None:
        ax.set_xticklabels(f'{t:.0f}' if t != 600 else r'$\infty$' for t in ax.get_xticks())
        #ax.set_xlabel(rf"Simulated ${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{process}}}}}$ [GeV]", usetex=True)
        ax.set_xlabel(rf"Simulated $\mathsf{{p}}_{{\text{{T}}}}^{{\text{{{process}}}}}$ \raisebox{{0.25ex}}{{[}}$\text{{GeV}}$\raisebox{{0.25ex}}{{]}}", usetex=True)
    else:
        ax_r.set_xticklabels(f'{t:.0f}' if t != 600 else r'$\infty$' for t in ax.get_xticks())
        #ax_r.set_xlabel(rf"Simulated ${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{process}}}}}$ [GeV]", usetex=True)
        ax_r.set_xlabel(rf"Simulated $\mathsf{{p}}_{{\text{{T}}}}^{{\text{{{process}}}}}$ \raisebox{{0.25ex}}{{[}}$\text{{GeV}}$\raisebox{{0.25ex}}{{]}}", usetex=True)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_yticks([10,100,1000,10000],minor=False)
    if doLog:
        ax.set_yscale('Log')
    scal_er = 4 if process == 'Z' else 4
    bott_er = 5 if process == 'Z' else 2
    #ax.set_ylim(2, 7000)
    ax.set_ylim(2, 8000)
    #ax.set_yticklabels(f'{t:.0f}' if t != 1 else '' for t in ax.get_yticks())
    if errorbox:
        handles, labels = ax.get_legend_handles_labels()
        from matplotlib import lines
        extra_patch = Patch(label='95\% CL interval',  fc='w', snap=True, alpha=0.0)
        seight_patch = Patch(label='68\% Expected',  fc=brazil_green, snap=True, edgecolor='k', linewidth=.3, alpha=1)
        nfive_patch = Patch(label='95\% Expected',  fc=brazil_yellow, snap=True, edgecolor='k', linewidth=.3, alpha=1)
        pred_patch = Patch(label='SM',  fc='#ea7dfd', snap=True, alpha=0.3)
        obs_hline = lines.Line2D([], [],  linestyle='-', color='k', 
                                 markersize=10, markeredgewidth=1.5, solid_capstyle='butt')
        exp_hline = lines.Line2D([], [],  linestyle='--', color='#046cf4', linewidth=0.5,
                                 markersize=10, markeredgewidth=1.5, solid_capstyle='butt')
        pred_hline = lines.Line2D([], [],  linestyle='-', color='#ea7dfd', 
                                  markersize=10, markeredgewidth=1.5, solid_capstyle='butt')
        #handles = [(asi_patch,handles[0]),(theo_patch,handles[0])]#[(asi_patch,handles[0])]
        handles = [obs_hline,exp_hline,seight_patch,nfive_patch,obs_hline,(pred_patch,pred_hline)]
        #labels  = ['Observed','Median expected','68% expected','95% expected']
        labels  = ['Obs.','Median exp.','68\% exp.','95\% exp.']
        leg_nf = ax.legend(handles,labels, #bbox_to_anchor=(1.00,1), 
                           fontsize=8, framealpha = 0, loc='upper right', handlelength=1.2, title='95\% CL upper limits', title_fontsize=8, ncol=2)
        #p_type_latex = { 'ttH' : r't$\mathregular{\bar{t}}$H', 'ttZ': r't$\mathregular{\bar{t}}$Z'}
        p_type_latex = { 'ttH' : r'$\mathsf{t\bar{t}H}$', 'ttZ': r'$\mathsf{t\bar{t}Z}$'}
        ax.legend([(pred_patch,pred_hline)],[f'SM {p_type_latex[p_type]}'], #bbox_to_anchor=(1.00,1), 
                  fontsize=8, framealpha = 0, loc='upper left', handlelength=1.2)
        ax.add_artist(leg_nf)
    else:
        plt.legend()

@save_pdf("exp_diffxsec_limits_lowptfrozen_ratio_final.pdf")
#@save_pdf("exp_diffxsec_limits_lowptfrozen_ratio.pdf")
def make_exp_xsec(p_xs, p_asymlim):
    bins = {'inc':np.array([200,600]),'stxs':np.array([200,300,450,600])}
    for p in p_xs:
        plot_exp_xsec(p_xs[p], p_asymlim, bins, p)

def get_p_xs():
    _out_dict = {}
    p_df = json.load(open(cfg.dataDir+'/process_norms/process_norms_ttbbw_run2.json','r'))
    # also get mu_rf, mu_r, mu_f, pdf, alphas
    scale = ['mu_f','mu_r','mu_rf']
    pdfas = ['pdfweight','alphas']
    # ttH scale uncertainty
    ttH_scale = {'Up': [0.89815028, 0.8827226 , 0.86699575, 0.84132501],
                 'Down':   [1.08194549, 1.10791701, 1.13462957, 1.18136289]}
    for p in ['ttZ','ttH']:
        o_df = {**{p:0},**{p+str(i):{'yield':0, 'theo_Up':0, 'theo_Down':0} for i in range(4)}}
        ttH_theo_by_year = {'2016':{}, '2017':{}, '2018': {}}
        for y in ['2016','2018','2017']:
            o_df[p] = o_df[p]+p_df[y][p]['stxs_yield']
            for i in range(4):
                o_df[p+str(i)]['yield'] = o_df[p+str(i)]['yield']+p_df[y][p+str(i)]['stxs_yield']
                for ud_str in ['Up','Down']:
                    if p == 'ttH':
                        scale_unc = abs(1-ttH_scale[ud_str][i])
                    else:
                        scale_unc = max([abs(1-1/p_df[y][p+str(i)][f'{sc}_{ud_str}']) for sc in scale])
                    pdfas_unc = np.sqrt((1-1/p_df[y][p+str(i)][f'pdfweight_{ud_str}'])**2 + (1-1/p_df[y][p+str(i)][f'alphas_{ud_str}'])**2)
                    theo_unc = np.sqrt((scale_unc)**2 + (pdfas_unc)**2)
                    if ud_str == 'Up':
                        theo_unc = 1+theo_unc
                    else:
                        theo_unc = 1-theo_unc
                    # because 2017 ttH theo is not supported 
                    if p == 'ttH':
                        ttH_theo_by_year[y][ud_str] = theo_unc
                    if y == '2017' and p == 'ttH':
                        theo_unc = np.mean([ttH_theo_by_year['2016'][ud_str], ttH_theo_by_year['2018'][ud_str]])
                    o_df[p+str(i)][f'theo_{ud_str}'] = o_df[p+str(i)][f'theo_{ud_str}']+p_df[y][p+str(i)]['stxs_yield']*theo_unc

            #
        #
        xs_df = {**{f'{p}':signal_xsec[p]*1000}, **{f'{p}{i}': {k:signal_xsec[p]*o_df[p+str(i)][k]/o_df[p]*1000 for k in o_df[p+str(i)]} for i in range(4)}}
        _out_dict[p] = xs_df
    return _out_dict

def parse_asym_limits():
    _out_dict = {}
    for k in asym_limits: # each signal strength modifier
        lims = asym_limits[k].split('\n')[1:]
        _out_dict[re.search(r'tt(Z|H)\d*',k).group()] = {
            'obs':float(lims[0].split('<')[-1].strip()),
            'exp':np.array([float(l.split('<')[-1].strip()) for l in lims[1:]])
        }
    return _out_dict

def print_latex_format(y_obs, y_exp, y_theo, sigp): # formats and prints limits in latex format
    y_exp_psig = y_exp[:,-2]-y_exp[:,2]
    y_exp_msig = y_exp[:,2]-y_exp[:,1]
    pt_dict = {0:"(200, 300]", 1:"(300, 450]", 2:r"(450, $\infty$)"}
    print("")
    for i in range(len(y_exp)):
        out_string = ""
        if i == 0:
            out_string += rf"\{sigp} "
        out_string += rf"& {pt_dict[i]:15} & ${y_obs[i]:.3f}$ (${{{y_exp[i,2]:.3f}}}_{{-{y_exp_msig[i]:.3f}}}^{{+{y_exp_psig[i]:.3f}}}$)"
        out_string +=rf" & ${y_obs[i]/y_theo[i]:.3f}$ (${{{y_exp[i,2]/y_theo[i]:.3f}}}_{{-{y_exp_msig[i]/y_theo[i]:.3f}}}^{{+{y_exp_psig[i]/y_theo[i]:.3f}}}$) \\" 
        if i == len(y_exp) and sigp == 'ttZ':
            out_string += r"[\cmsTabSkip]"
        print(out_string)
    
def main():
    import_mpl_settings(1, length=1.25)
    p_xs = get_p_xs()
    p_asymlims = parse_asym_limits()
    #
    make_exp_xsec(p_xs, p_asymlims)

if __name__ == '__main__':
    main()
