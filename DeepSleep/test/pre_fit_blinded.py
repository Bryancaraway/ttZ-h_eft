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
from lib.fun_library import save_pdf, getLaLabel
from qcDatacard import tbins_map
from post_fit import PostFit
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            




@save_pdf('blinded_prefit.pdf')
def main():
    froo = f'fitDiagnostics_blind_run2.root'
    prefit= PreFit(froo)
    prefit.makeplots(doPull=False)

class PreFit(PostFit):
    '''
    steal some of the functionality 
    from already existing class
    '''
    def __init__(self, fitroo):
        super().__init__(fitroo)


    def do_stack(self, d,ax,ch):
        ycumm = None
        # stack portion
        #ordered_list = re.findall(rf'tt[H,Z]\d', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_2b','tt_bb','TTBar']
        ordered_list = re.findall(rf'tt[H,Z]\d', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_B','TTBar']
        #colors =  plt.cm.gist_rainbow(np.linspace(0,1,len(ordered_list)))
        colors =  plt.cm.tab20c(np.linspace(0,1,20))[:8]
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

    def make_stackhist(self, ch='Zhpt2'):
        ax = self.axs # # 0,0 lt # 0,1 rt
        d = self.hists['prefit'][ch]
        ## stack plot here
        self.do_stack(d, ax, ch)
        #
        self.endaxs(ax,'prefit', ch=ch)
        self.addlegend(opt='blind')
        #plt.show()

    def init_axes(self, opt='', channel=''):
        #fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig, axs = plt.subplots(1,1)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.0,
            wspace=0.3
        )
        self.fig = fig
        lumi = cfg.Lumi[re.search(r'201\d',channel).group()]
        self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Simulation$", fontsize = 10)
        self.fig.text(0.70,0.89, f'{lumi:.1f}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
        self.axs = axs

    def endaxs(self,ax,p,ch='Zhpt1'):
        #axt
        ax.set_xlim(0,self.edges[ch][-1])
        ax.set_ylim(0.)

        #if axt.get_ylim()[1] > 1000:
        ax.set_yscale('log')
        ax.set_ylim(ymin=.1, ymax=ax.get_ylim()[1] * 10)
        #else:
        #     axt.yaxis.set_minor_locator(AutoMinorLocator())

        ax.set_ylabel('Events / bin', fontsize=8)
        ax.set_title(rf'${p}$', y=1.0, pad=-14, fontsize=8)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

        ax.tick_params(which='both', direction='in', top=True, right=True)

if __name__ == '__main__':
    main()
