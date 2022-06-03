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


#fit_vars = ['Zh_bbvLscore', 'n_b_inZh','outZH_b1_q_mindr',]
k_format = {
    'Zh_bbvLscore' : r'$\mathsf{Z/H\;b\bar{b}\;score}$',
    'n_b_inZh': r'$N\mathsf{(b_{in})}$',
    'outZH_b1_q_mindr': r'$\mathsf{\Delta} R \mathsf{(b,q)}$',
    'n_ak8jets': r'$N\mathsf{(AK8\;jets)}$',
    'ak4_worstb_inZH': r'$\mathsf{min.\;b_{in}\;score}$',
    'inZhb_outZhb_dr': r'$\mathsf{\Delta} R \mathsf{(b_{in},b_{out})}$',
}
b_str = {
    'Zh_bbvLscore': '/ 0.25', 
    'n_b_inZh': '', 
    'outZH_b1_q_mindr': '/ 0.4',
    'n_ak8jets':'',
    'ak4_worstb_inZH' :'/ 0.1',
    'inZhb_outZhb_dr':'bin',
}
tick_dict = {
    'n_b_inZh': [0,1,2],
    'n_ak8jets': [1,2,3],
    'outZH_b1_q_mindr': [0.8,1.6,2.4,3.2],
    'inZhb_outZhb_dr': [0.8,1.6,2.4,3.2],
}
#pre_or_post = 'postfit'
pre_or_post = 'prefit'

@save_pdf('nn_pas_validation.pdf')
def main():
    for v in k_format:
        print(v)
        froo = f'fitDiagnostics_{v}_NNcuts_run2.root'
        pfp= PlotForPasNN(froo,v,tbins_map[v])
        pfp.makeplots()

class PlotForPasNN(PostFit):
    '''
    steal some of the functionality 
    from already existing class
    '''
    def __init__(self, fitroo, kinem, t_labels):
        super().__init__(fitroo, kinem=kinem, t_labels=t_labels)

    def makeplots(self):
        ch_list = list(self.hists[pre_or_post].keys())
        hist_list = list(self.hists[pre_or_post][ch_list[0]].keys())
        for hist_n in hist_list:
            self.hists[pre_or_post][ch_list[-1]][hist_n]['values'] = np.sum([self.hists[pre_or_post][ch][hist_n]['values'] for ch in ch_list], axis=0)
            if hist_n != 'data':
                self.hists[pre_or_post][ch_list[-1]][hist_n]['err'] = np.sqrt(np.sum([np.power(self.hists[pre_or_post][ch][hist_n]['err'],2) for ch in ch_list], axis=0))
            else:
                self.hists[pre_or_post][ch_list[-1]][hist_n]['errup'] = np.sqrt(np.sum([np.power(self.hists[pre_or_post][ch][hist_n]['errup'],2) for ch in ch_list], axis=0))
                self.hists[pre_or_post][ch_list[-1]][hist_n]['errdw'] = np.sqrt(np.sum([np.power(self.hists[pre_or_post][ch][hist_n]['errdw'],2) for ch in ch_list], axis=0))
        self.init_axes(year='2018') # dummy year
        self.make_stackhist(ch_list[-1], year='2018') # dummy channel and year


    #def makeplots(self): # year-by-year plots
    #    for ch in self.hists[pre_or_post]:
    #        y = re.search(r'201\d',ch).group()
    #        self.init_axes(year=y)
    #        self.make_stackhist(ch, year=y)

    def do_stack(self, d,ax,ch):
        ycumm = None
        # stack portion
        ordered_list = ['ttZ','ttH','VJets','ttX','single_t','tt_B','TTBar']
        colors =  plt.cm.tab10(np.linspace(0,1,10))[0:2]
        colors = np.append(colors, plt.cm.gist_rainbow(np.linspace(0,1,6)), axis=0)
        #
        for j,k in enumerate(ordered_list):
            if k not in d: continue
            y = np.append(d[k]['values'],0)
            if ycumm is not None:
                ycumm += y
            else:
                ycumm = y 
            c = colors[j]
            label = getLaLabel(k)[0]
            ax.fill_between(self.edges[ch],ycumm,ycumm-y,step="post", 
                            linewidth=0, color=c, label=label)
            #
        #
        self.make_error_boxes(ax, (self.edges[ch][1:]+self.edges[ch][:-1])/2, d['total']['values'],
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err'], label='Stat.+syst.')
        
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
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err']/d['total']['values'], label='Stat.+syst.')
        ax.set_ylim(0.25,1.75)
        ax.set_yticks([.5,1,1.5])
        ax.tick_params(which='both', direction='in', bottom=True, left=True, right=True, top=True) 
        ax.grid(True)
        ax.set_xlim(self.edges[ch][0],self.edges[ch][-1])        

    def make_stackhist(self, ch, year):
        ax, rat_ax = self.ax, self.rat_ax
        d = self.hists[pre_or_post][ch]
        ## stack plot here
        self.do_stack(d, ax, ch)
        self.do_data(d, ax, ch)
        self.do_ratio(d, rat_ax, ch)
        self.endaxs(ax, pre_or_post, ch=ch)
        #
        b_w = (self.edges[ch][1:]-self.edges[ch][:-1])/2
        ax.set_ylabel(f'Events {b_str[self.kinem]}')#, fontsize=8)
        ax.set_yscale('log')
        ax.set_ylim(ymin=1.01, ymax=max(d['total']['values'])*600)#ymax=ax.get_ylim()[1] * 500)
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, self.fig.dpi_scale_trans)
        ext_args = {'usetex':True, 'transform':trans, 'fontsize':8.5}
        #ax.text(0.03, .90, rf"{year} $\textit{{{'P'+pre_or_post[1:]}}}$" , **ext_args)
        ax.text(0.03, .90, rf"$\textit{{{'P'+pre_or_post[1:]}}}$" , **ext_args)
        #
        rat_ax.set_xlabel(f'{k_format[self.kinem]}', labelpad=0.0)
        rat_ax.set_ylabel('Data / MC')
        self.addlegend(opt='blind')
        #plt.tight_layout(h_pad=0)
        #plt.show()

    def init_axes(self, opt='', year=''):
        fig, (ax,rat_ax) = plt.subplots(2,1, sharey='row', sharex='col', gridspec_kw={'height_ratios':[4,1]})
        fig.subplots_adjust(
            top=0.95,
            #top=0.88,
            bottom=0.12,
            left=0.15,
            right=0.95,
            hspace=0.0,
            wspace=0.0
        )
        self.fig = fig
        lumi = cfg.Lumi[re.search(r'201\d',year).group()]
        #CMSlabel(fig, ax, lumi=lumi, altax=ax, opt='')
        CMSlabel(fig, ax, altax=ax, opt='')
        self.ax = ax
        self.rat_ax = rat_ax

    def endaxs(self,ax,p,ch='Zhpt1'):
        #axt
        ax.set_xlim(self.edges[ch][0],self.edges[ch][-1])
        ax.tick_params(which='both', direction='in', bottom=True, left=True, top=True, right=True)
        if 'dr' in self.kinem:
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        else:
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        if self.kinem in tick_dict:
            ax.set_xticks(tick_dict[self.kinem])

    def addlegend(self, opt=''):
        ax = self.ax
        handles, labels = ax.get_legend_handles_labels()
        hatch_patch = Patch(hatch=10*'X', label='Stat.+syst.',  fc='w', alpha=0.99)
        #handles = handles[0:2]+handles[-3:-1]+handles[2:-3] + [handles[-1]] + [hatch_patch]
        #labels  = labels[0:2]+labels[-3:-1]+labels[2:-3] + [labels[-1]] + ['Stat+sys.']
        handles = handles[:-1][::-1] + [handles[-1]] +[hatch_patch]
        labels = labels[:-1][::-1] + [labels[-1]] +['Stat.+syst.']
        ax.legend(handles,labels, bbox_to_anchor=(1.02, 1.03), 
                  ncol = 2, #columnspacing=8.0,
                  fontsize=8.5, framealpha = 0, loc='upper right')

if __name__ == '__main__':
    import_mpl_settings(1, length=1.25)
    main()



