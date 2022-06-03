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


#pre_or_post = 'postfit'
pre_or_post = 'prefit'

wc_latex = {
    'ctp'  : r'$\mathsf{c_{t \varphi}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cpQM' : r'$\mathsf{c^{-}_{\varphi Q}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cpQ3' : r'$\mathsf{c^3_{\varphi Q}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cpt'  : r'$\mathsf{c_{\varphi t}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cptb' : r'$\mathsf{c_{\varphi t b}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'ctW'  : r'$\mathsf{c_{tW}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cbW'  : r'$\mathsf{c_{bW}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'ctZ'  : r'$\mathsf{c_{tZ}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
}

#@save_pdf('pas_eftfit_plots_unblind_test.pdf')
@save_pdf('pas_eftfit_plots_unblind_prefit_ttbb.pdf')
#@save_pdf('pas_eftfit_plots_unblind_postfit_ttbb.pdf')
#@save_pdf('pas_eftfit_plots_unblind_fixed_testformat.pdf')
def main():
    froo = 'fitDiagnostics_unblind_SM_run2.root'# unblind fixed to SM
    wc_ranges ={
        #'ctW'  :[-1.03, 0.94]  ,
        #'ctZ'  :[-0.99, 1.02]  ,
        #'ctp'  :[0.13, 30.18],
        #'cpQM' :[-4.77, 5.63],
        ##'ctG'  :[-0.50, 0.48]  ,
        #'cbW'  :[-4.55, 4.65]  ,
        #'cpQ3' :[-3.86, 2.88]  ,
        #'cptb' :[-9.44, 10.15],
        #'cpt'  :[-8.50, 5.39]
        'ctW'  :[-1.05, 0.96, -0.05]  ,
        'ctZ'  :[-1.05, 1.11, 0.03]  ,
        'ctp'  :[0.25, 30.04, 15.25],
        'cpQM' :[-6.55, 8.72, 0.86],
        #'ctG'  :[-0.50, 0.48]  ,
        'cbW'  :[-5.05, 5.08, -2.00]  ,
        'cpQ3' :[-4.13, 3.04, -0.50]  ,
        'cptb' :[-9.88, 10.75, 0.58],
        'cpt'  :[-12.02, 6.32, -2.09]
    }
    for wc in wc_ranges:
        pfp= EFTFitForPas(froo, wc, wc_ranges[wc])
        pfp.makeplots()

class EFTFitForPas(PostFit):
    '''
    steal some of the functionality 
    from already existing class
    '''
    def __init__(self, fitroo, wc, wc_range):
        self.wc = wc
        self.wc_range = wc_range
        super().__init__(fitroo)


    def makeplots(self):
        #self.do1Pull = do1Pull
        ch_list = [ch for ch in self.hists[pre_or_post]] # 'prefit is just a placeholder to loop over pt'
        ch_groups = {y: re.findall(rf'\w*{y}\w*' ,' '.join(ch_list)) for y in cfg.Years}
        #self.fig.suptitle(ch, fontsize=10) # maybe get rid of 
        for y in cfg.Years:
            self.eft_ycumm = {}
            self.init_axes(year=y)
            self.make_stackhist(ch_groups[y])

    def do_stack(self, d,ax,ch):
        ycumm = None
        ordered_list = ['ttZ','ttH','VJets','ttX','single_t','tt_B','TTBar']
        d['ttZ'] = {'values' : np.nansum([d[f'ttZ{i}']['values'] for i in range(4)], axis=0)}
        d['ttH'] = {'values' : np.nansum([d[f'ttH{i}']['values'] for i in range(4)], axis=0)}
        #
        for j,k in enumerate(ordered_list):
            #if 'data' in k or 'total' in k: continue
            if k not in d: continue
            y = np.append(d[k]['values'],0)
            if ycumm is not None:
                ycumm += y
            else:
                ycumm = y 
        self.ycumm = ycumm
        
    def do_data(self, d, ax, ch):
        ax.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=d['data']['values'],
                    xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,yerr=[d['data']['errdw'],d['data']['errup']], 
                    fmt='.', label='Data', color='k')
    def do_ratio(self, d, ax, ch):
        for wc_v in self.eft_ycumm:
            y       =  self.eft_ycumm[wc_v] / self.ycumm
            ax.step(x=self.edges[ch], y=y,
                    where='post',
                    label=f'{self.wc} = {self.wc_range[wc_v]}',)
                    #color='tab:blue')
        self.make_error_boxes(ax, (self.edges[ch][1:]+self.edges[ch][:-1])/2, np.ones_like(d['total']['values']),
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err']/d['total']['values'], label='stat+sys')
        ax.set_ylim(0.5,2.0)
        ax.set_yticks([1,1.5])
        if 'Zhpt1' in ch:
            ax.tick_params(which='both', direction='in', bottom=True, left=True, right=False, top=False) 
        else:
            ax.tick_params(which='both', direction='in', bottom=True, left=False, right=False, top=False) 
        ax.grid(True)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(0,self.edges[ch][-1])

    def do_vlines(self, ax, i_, ch, isratio=False):
        m_factor = 3 if i_ == 0 else 4
        m_lines = [ self.edges[ch][i] for i in range(m_factor,len(self.edges[ch]), m_factor)]
        for l in m_lines:
            ymax = 0.60 if i_ != 0  else 0.8
            ymax = 1.0 if isratio else ymax
            ax.axvline(l, ymin=-.05, ymax=ymax,
                       color='k', linewidth='0.5', linestyle='--', dashes=(4,8), snap=True)
        

    def make_stackhist(self, ch_list):
        axs = self.axs # # 0,0 lt # 0,1 rt
        #rat_axs = self.rat_axs
        for i, ch in enumerate(ch_list):
            d = self.hists[pre_or_post][ch]
            ## stack plot here
            self.do_stack(d, axs[i], ch)
            #self.do_data(d, axs[i], ch)
            self.compute_eft_impact(d, ch, 1)
            self.compute_eft_impact(d, ch, 2)
            self.do_ratio(d, axs[i], ch)
            self.do_vlines(axs[i], i, ch)
            #self.do_vlines(rat_axs[i], i, ch, isratio=True)
            #
            self.endaxs(axs[i],pre_or_post, ch=ch)
        #
        axs[0].set_ylabel('(SM+EFT)/SM', fontsize=14)
        #axs[0].set_yscale('log')
        #axs[0].set_ylim(ymin=.05, ymax=axs[0].get_ylim()[1] * 10)
        #axs[0].set_ylim(ymin=.05, ymax=axs[0].get_ylim()[1] * 2)
        axs[-1].set_xlabel('Analysis bins',fontsize=14)
        #rat_axs[0].set_ylabel('Data / MC', fontsize=14)
        self.nn_text_bins(axs[0])
        self.addlegend(opt='blind')
        #plt.tight_layout(h_pad=0,w_pad=0)
        #plt.show()


    def nn_text_bins(self, ax):
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, self.fig.dpi_scale_trans)
        ext_args = {'usetex':True, 'transform':trans}
        ax.text(0.05, .92, rf"$\textit{{{'P'+pre_or_post[1:]}}}$" , **ext_args)
        ax.text(0.05, .85, r'NN${}_{\mathrm{bin}}$' , **ext_args)
        ax.text(0.04, .80, '1', **ext_args)
        ax.text(0.215, .80, '2', **ext_args)
        ax.text(0.39, .80, '3', **ext_args)

        ax.text(0.565, .80, '4', **ext_args)
        ax.text(0.73, .80, '5', **ext_args)
        ax.text(0.895, .80, '6', **ext_args)


    def init_axes(self, opt='', year=''):
        #fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        #fig, axs = plt.subplots(1,3, sharey=True)
        fig, axs = plt.subplots(1,3, sharey='row', sharex='col', )
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
        #CMSlabel(fig, axs[0], lumi=lumi, altax=axs[-1], fontsize=14, opt='')
        CMSlabel(fig, axs[0], lumi=lumi, altax=axs[-1], fontsize=14)
        self.axs = axs
        #self.rat_axs = rat_axs

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
            ax.tick_params(which='both', direction='in', bottom=True, left=True, top=False, right=False)
        else:
            ax.tick_params(which='both', direction='in', bottom=True, left=True, top=False, right=False)
        ax.tick_params(which='both', direction='in')
        channel = re.search(r'Zhpt\d',ch).group()
        ax.set_xticks(tick_dict[channel][0])
        ax.set_xticklabels(tick_dict[channel][1])
        #ax.set_title(f"{ch_dict[channel]}", loc='right', y=1.0, pad=-15, fontsize=10, usetex=True)
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, self.fig.dpi_scale_trans)
        ext_args = {'usetex':True, 'transform':trans}
        ax.text( .95, .92, f"{ch_dict[channel]}", ha='right', **ext_args)
        #ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    def addlegend(self, opt=''):
        #ax = self.axs[0,1]
        ax = self.axs[-2]
        handles, labels = ax.get_legend_handles_labels()
        print(labels)
        hatch_patch = Patch(hatch=10*'X', label='stat+sys',  fc='w', alpha=0.99)
        handles = handles + [hatch_patch]
        wc_labels = [wc_latex[label.split('=')[0].strip()] + ' = ' +label.split('=')[1].strip() for label in labels]
        labels  = wc_labels + ['stat+sys.']
        ax.legend(handles,labels, bbox_to_anchor=(.3205,.740), 
                  ncol = 1, #columnspacing=4.0,
                  fontsize=10, framealpha = 0, loc='lower left')

    def compute_eft_impact(self, d, ch, index):
        fit_f   = pd.read_pickle(sys.path[1]+"/Higgs-Combine-Tool/EFT_Parameterization_v7.npy")
        #('tt_B', 'y2016_Zhpt1_0')
        get_eft_w = (lambda s_,v_ : np.array([fit_f[(s_,ch+'_'+str(i_))][(self.wc, self.wc)] * self.wc_range[v_] * self.wc_range[v_] +\
                                              fit_f[(s_,ch+'_'+str(i_))][('sm', self.wc)] * self.wc_range[v_] +\
                                              fit_f[(s_,ch+'_'+str(i_))][('sm', 'sm')]
                                              for i_ in range(len(self.edges[ch])-1)]))
        samples = ['ttZ','ttH','VJets','ttX','single_t','tt_B','TTBar']
        eft_ycumm = None
        #print(self.ycumm)
        for s in samples:
            #if s in ['tt_B','ttZ','ttH']:
            #if s in ['ttZ','ttH']:
            if s in ['tt_B']:
                y = np.append(d[s]['values']*get_eft_w(s,index),0)
            else:
                y = np.append(d[s]['values'],0)
            if eft_ycumm is not None:
                eft_ycumm += y
            else:
                eft_ycumm = y
        self.eft_ycumm[index] = eft_ycumm
    

if __name__ == '__main__':
    import_mpl_settings(1, length=2, width=3, disable_sansmath=True)
    main()
