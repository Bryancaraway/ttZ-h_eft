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
import re
import seaborn as sns
import config.ana_cff as cfg
#from lib.fun_library import save_pdf, getLaLabel
from lib.fun_library import t2Run, save_pdf, getZhbbBaseCuts, getZhbbWeight, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel
from qcDatacard import tbins_map
from post_fit import PostFit
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import transforms
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
#rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            


pre_or_post = 'postfit'
#pre_or_post = 'prefit'


#@save_pdf('pas_template_plots_unblind_SMprefit.pdf')
#@save_pdf('pas_template_plots_unblind_SM.pdf')
@save_pdf('pas_template_plots_unblind_SM_final.pdf')
def main():
    #froo = f'fitDiagnostics_inc_run2.root'
    #froo = f'fitDiagnostics_inctest_run2.root'
    ###
    froo = 'fitDiagnostics_unblind_SM_run2.root'# unblind fixed to SM
    #froo = 'fitDiagnostics_unblind_mu0_run2.root'# unblind constrained > 0
    pfp= PlotForPas(froo)
    pfp.makeplots()

class PlotForPas(PostFit):
    '''
    steal some of the functionality 
    from already existing class
    '''
    def __init__(self, fitroo):
        super().__init__(fitroo)


    def makeplots(self):
        #self.do1Pull = do1Pull
        ch_list = [ch for ch in self.hists[pre_or_post]] # 'prefit is just a placeholder to loop over pt'
        ch_groups = {y: re.findall(rf'\w*{y}\w*' ,' '.join(ch_list)) for y in cfg.Years}
        #self.fig.suptitle(ch, fontsize=10) # maybe get rid of 
        for y in cfg.Years:
            self.init_axes(year=y)
            self.make_stackhist(ch_groups[y], year=y)

    def do_stack(self, d,ax,ch):
        ycumm = None
        # stack portion
        #ordered_list = re.findall(rf'tt[H,Z]\d', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_B','TTBar']
        ordered_list = ['ttZ','ttH','VJets','ttX','single_t','tt_B','TTBar']
        d['ttZ'] = {'values' : np.nansum([d[f'ttZ{i}']['values'] for i in range(4)], axis=0)}
        d['ttH'] = {'values' : np.nansum([d[f'ttH{i}']['values'] for i in range(4)], axis=0)}
        #
        #colors =  plt.cm.tab20c(np.linspace(0,1,20))[:8]
        colors =  plt.cm.tab10(np.linspace(0,1,10))[0:2]
        colors = np.append(colors, plt.cm.gist_rainbow(np.linspace(0,1,6)), axis=0)
        #
        for j,k in enumerate(ordered_list):
            #if 'data' in k or 'total' in k: continue
            if k not in d: continue
            y = np.append(d[k]['values'],0)
            if ycumm is not None:
                ycumm += y
            else:
                ycumm = y 
            c = colors[j]
            #c = colors[j + (len(colors)//2)*(j % 2) - 1*((j+1)%2)]
            label = getLaLabel(k)[0]
            ax.fill_between(self.edges[ch],ycumm,ycumm-y,step="post", 
                            linewidth=0, color=c, label=label)
            #
        #
        self.make_error_boxes(ax, (self.edges[ch][1:]+self.edges[ch][:-1])/2, d['total']['values'],
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err'], label='Stat+sys')
        
    def do_data(self, d, ax, ch):
        ax.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=d['data']['values'],
                    xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,yerr=[d['data']['errdw'],d['data']['errup']], 
                    fmt='.', label='Data', color='k')
    def do_ratio(self, d, ax, ch):
        y       =  d['data']['values']/d['total']['values']
        yerrdw  =  d['data']['errdw'] /d['total']['values']
        yerrup  =  d['data']['errup'] /d['total']['values']
        ax.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=y,
                     xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,yerr=[yerrdw,yerrup], 
                     fmt='.', label='Data', color='k')
        self.make_error_boxes(ax, (self.edges[ch][1:]+self.edges[ch][:-1])/2, np.ones_like(d['total']['values']),
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err']/d['total']['values'], label='Stat+sys')
        ax.set_ylim(0,2)
        ax.set_yticks([.5,1,1.5])
        if 'Zhpt1' in ch:
            ax.tick_params(which='both', direction='in', bottom=True, left=True, right=False, top=False, labelsize=12) 
        else:
            ax.tick_params(which='both', direction='in', bottom=True, left=False, right=False, top=False, labelsize=12) 
        ax.grid(True)
        #ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(0,self.edges[ch][-1])

    def do_vlines(self, ax, i_, ch, isratio=False):
        m_factor = 3 if i_ == 0 else 4
        m_lines = [ self.edges[ch][i] for i in range(m_factor,len(self.edges[ch]), m_factor)]
        for l in m_lines:
            ymax = 0.55 if i_ != 0  else 0.8
            ymax = 1.0 if isratio else ymax
            ax.axvline(l, ymin=-.05, ymax=ymax,
                       color='k', linewidth='0.5', linestyle='--', dashes=(4,8), snap=True)
        

    def make_stackhist(self, ch_list, year):
        axs = self.axs # # 0,0 lt # 0,1 rt
        rat_axs = self.rat_axs
        for i, ch in enumerate(ch_list):
            d = self.hists[pre_or_post][ch]
            ## stack plot here
            self.do_stack(d, axs[i], ch)
            self.do_data(d, axs[i], ch)
            self.do_ratio(d, rat_axs[i], ch)
            self.do_vlines(axs[i], i, ch)
            self.do_vlines(rat_axs[i], i, ch, isratio=True)
            #
            self.endaxs(axs[i],pre_or_post, ch=ch)
        #
        axs[0].set_ylabel('Events / bin', fontsize=14)
        axs[0].set_yscale('log')
        axs[0].set_ylim(ymin=.05, ymax=axs[0].get_ylim()[1] * 10)
        rat_axs[-1].set_xlabel(f'{year} analysis bins',fontsize=14)
        rat_axs[0].set_ylabel('Data / MC', fontsize=14)
        self.nn_text_bins(axs[0], year)
        self.addlegend(opt='blind')
        #plt.tight_layout(h_pad=0,w_pad=0)
        #plt.show()


    def nn_text_bins(self, ax, year):
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, self.fig.dpi_scale_trans)
        ext_args = {'usetex':True, 'transform':trans, 'fontsize':12}
        ax.text(0.04, .92, rf"$\textit{{{'P'+pre_or_post[1:]}}}$" , **ext_args)
        ax.text(0.04, .83, r'DNN${}_{\mathrm{bin}}$' , **ext_args)
        ax.text(0.07, .78, '1', **ext_args)
        ax.text(0.23, .78, '2', **ext_args)
        ax.text(0.40, .78, '3', **ext_args)

        ax.text(0.57, .78, '4', **ext_args)
        ax.text(0.73, .78, '5', **ext_args)
        ax.text(0.90, .78, '6', **ext_args)


    def init_axes(self, opt='', year=''):
        #fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        #fig, axs = plt.subplots(1,3, sharey=True)
        fig, (axs,rat_axs) = plt.subplots(2,3, sharey='row', sharex='col', gridspec_kw={'height_ratios':[4,1]})
        fig.subplots_adjust(
            top=0.95,
            bottom=0.11,
            left=0.11,
            right=0.95,
            hspace=0.0,
            wspace=0.0
        )
        self.fig = fig
        lumi = cfg.Lumi[re.search(r'201\d',year).group()]
        #self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Simulation$", fontsize = 10)
        #self.fig.text(0.70,0.89, f'{lumi:.1f}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
        CMSlabel(fig, axs[0], lumi=lumi, altax=axs[-1], fontsize=16, opt='')
        #CMSlabel(fig, axs[0], lumi=lumi, altax=axs[-1], fontsize=14)
        self.axs = axs
        self.rat_axs = rat_axs

    def endaxs(self,ax,p,ch='Zhpt1'):
        #axt
        ax.set_xlim(0,self.edges[ch][-1])
        #ax.set_ylim(0.)
        ch_dict = {
            #'Zhpt1': r'$200 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 300$ GeV ',
            #'Zhpt2': r'$300 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 450$ GeV ',
            #'Zhpt3': r'${p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} > 450$ GeV '
            'Zhpt1': r'$\mathsf{200} < \mathsf{p}_{\text{T}}^{\mathrm{Z/H\;cand.}} < \mathsf{300}$ GeV ',
            'Zhpt2': r'$\mathsf{300} < \mathsf{p}_{\text{T}}^{\mathrm{Z/H\;cand.}} < \mathsf{450}$ GeV ',
            'Zhpt3': r'$\mathsf{p}_{\text{T}}^{\text{Z/H\;cand.}} > \mathsf{450}$ GeV '
        }
        tick_dict = { 
            'Zhpt1': [np.arange(0,19,3),np.arange(0,19,3)],
            'Zhpt2': [np.arange(4,25,4),np.arange(22,43,4)],
            'Zhpt3': [np.arange(4,25,4),np.arange(46,67,4)],
        }
        #if axt.get_ylim()[1] > 1000:
        #else:
        #     axt.yaxis.set_minor_locator(AutoMinorLocator())
        if 'Zhpt1' in ch:
            ax.tick_params(which='both', direction='in', bottom=True, left=True, top=False, right=False, labelsize=12)
        else:
            ax.tick_params(which='both', direction='in', bottom=True, left=True, top=False, right=False, labelsize=12)
        ax.tick_params(which='both', direction='in')
        channel = re.search(r'Zhpt\d',ch).group()
        ax.set_xticks(tick_dict[channel][0])
        ax.set_xticklabels(tick_dict[channel][1])
        #ax.set_title(f"{ch_dict[channel]}", loc='right', y=1.0, pad=-15, fontsize=10, usetex=True)
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, self.fig.dpi_scale_trans)
        ext_args = {'usetex':True, 'transform':trans, 'fontsize':12}
        ax.text( .97, .92, f"{ch_dict[channel]}", ha='right', **ext_args)
        #ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    def addlegend(self, opt=''):
        #ax = self.axs[0,1]
        ax1 = self.axs[1] # second pane
        ax2 = self.axs[-1] # third/last pane
        handles, labels = ax2.get_legend_handles_labels()
        print(labels)
        hatch_patch = Patch(hatch=10*'X', label='Stat.+sys',  fc='w', alpha=0.99)
        #handles = handles[0:2]+handles[-3:-1]+handles[2:-3] + [handles[-1]] + [hatch_patch]
        #labels  = labels[0:2]+labels[-3:-1]+labels[2:-3] + [labels[-1]] + ['Stat+sys.']
        handles = handles[:-1][::-1] + [handles[-1]] +[hatch_patch]
        labels = labels[:-1][::-1] + [labels[-1]] +['Stat.+sys.']
        ax1.legend(handles[:4],labels[:4], bbox_to_anchor=(1.05,.930), 
                   #ncol = 2,
                  fontsize=12, framealpha = 0, loc='upper right')
        ax2.legend(handles[4:], labels[4:], bbox_to_anchor=(.35,.930), 
                   #ncol = 2,
                   fontsize=12, framealpha = 0, loc='upper left')

if __name__ == '__main__':
    import_mpl_settings(1, length=2, width=3)
    main()

