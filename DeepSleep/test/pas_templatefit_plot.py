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
import uproot
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


@save_pdf('pas_template_plots.pdf')
def main():
    froo = f'fitDiagnostics_inc_run2.root'
    #froo = f'fitDiagnostics_nomerge_run2.root'
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
        ch_list = [ch for ch in self.hists['prefit']] # 'prefit is just a placeholder to loop over pt'
        ch_groups = {y: re.findall(rf'\w*{y}\w*' ,' '.join(ch_list)) for y in cfg.Years}
        #self.fig.suptitle(ch, fontsize=10) # maybe get rid of 
        for y in cfg.Years:
            self.init_axes(year=y)
            self.make_stackhist(ch_groups[y])

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
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err'], label='stat+sys')
        
    def do_vlines(self, ax, i_, ch):
        m_factor = 3 if i_ == 0 else 4
        m_lines = [ self.edges[ch][i] for i in range(m_factor,len(self.edges[ch]), m_factor)]
        for l in m_lines:
            ax.axvline(l, ymin=-.05, ymax=(0.7 if i_ != 0 else 0.8), 
                       color='k', linewidth='0.5', linestyle='--', dashes=(4,8), snap=True)
        

    def make_stackhist(self, ch_list):
        axs = self.axs # # 0,0 lt # 0,1 rt
        for i, ch in enumerate(ch_list):
            d = self.hists['prefit'][ch]
            ## stack plot here
            self.do_stack(d, axs[i], ch)
            self.do_vlines(axs[i], i, ch)
            #
            self.endaxs(axs[i],'prefit', ch=ch)
        #
        axs[0].set_ylabel('Events / bin', fontsize=10)
        axs[0].set_yscale('log')
        axs[0].set_ylim(ymin=.1, ymax=axs[0].get_ylim()[1] * 10)
        axs[-1].set_xlabel('Analysis bins',fontsize=10)
        self.nn_text_bins(axs[0])
        self.addlegend(opt='blind')
        #plt.tight_layout()
        #plt.show()


    def nn_text_bins(self, ax):
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, self.fig.dpi_scale_trans)
        ext_args = {'usetex':True, 'transform':trans}
        ax.text(0.05, .95, r'$\textit{Prefit}$' , **ext_args)
        ax.text(0.05, .83, r'NN${}_{\mathrm{bin}}$' , **ext_args)
        ax.text(0.04, .80, '1', **ext_args)
        ax.text(0.215, .80, '2', **ext_args)
        ax.text(0.39, .80, '3', **ext_args)

        ax.text(0.565, .80, '4', **ext_args)
        ax.text(0.73, .80, '5', **ext_args)
        ax.text(0.895, .80, '6', **ext_args)


    def init_axes(self, opt='', year=''):
        #fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig, axs = plt.subplots(1,3, sharey=True)
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
        CMSlabel(fig, axs[0], lumi=lumi, altax=axs[-1])
        self.axs = axs

    def endaxs(self,ax,p,ch='Zhpt1'):
        #axt
        ax.set_xlim(0,self.edges[ch][-1])
        #ax.set_ylim(0.)
        ch_dict = {
            'Zhpt1': r'$200 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 300$',
            'Zhpt2': r'$300 < {p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} < 450$',
            'Zhpt3': r'${p}_{\mathrm{T}}^{\mathrm{Z/H\;cand.}} > 450$'
        }
        tick_dict = { 
            'Zhpt1': [np.arange(0,19,3),np.arange(0,19,3)],
            'Zhpt2': [np.arange(4,25,4),np.arange(22,43,4)],
            'Zhpt3': [np.arange(4,25,4),np.arange(46,67,4)],
        }
        #if axt.get_ylim()[1] > 1000:
        #else:
        #     axt.yaxis.set_minor_locator(AutoMinorLocator())
        #ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.tick_params(which='both', direction='in')
        channel = re.search(r'Zhpt\d',ch).group()
        ax.set_xticks(tick_dict[channel][0])
        ax.set_xticklabels(tick_dict[channel][1])
        ax.set_title(f"{ch_dict[channel]}", y=1.0, pad=-25, fontsize=10, usetex=True)
        #ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    def addlegend(self, opt=''):
        #ax = self.axs[0,1]
        ax = self.axs[-1]
        handles, labels = ax.get_legend_handles_labels()
        hatch_patch = Patch(hatch=10*'X', label='stat+sys',  fc='w')
        handles = handles + [hatch_patch]
        labels  = labels + ['stat+sys.']
        ax.legend(handles,labels, bbox_to_anchor=(-1.025,.685), 
                  ncol = 2, columnspacing=2.0,
                  fontsize=10, framealpha = 0, loc='lower left')

if __name__ == '__main__':
    import_mpl_settings(2)
    main()
