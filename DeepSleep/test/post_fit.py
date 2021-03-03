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
np.random.seed(1)
import re
import uproot
from scipy.stats import norm
import config.ana_cff as cfg
from lib.fun_library import save_pdf, getLaLabel
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            


def main():
    # fitDiagnostics_test_2016.root, fitDiagnostics_test_2017.root, fitDiagnostics_test_2018.root, sfitDiagnostics_test_run2.root
    #PostFit("fitDiagnostics_test_2016.root").makeplots()
    #f_list = [f"fitDiagnostics_blind_{i}.root" for i in ['2016','2017','2018','run2']] # test
    f_list = ['fitDiagnostics_partblind_run2.root']
    #f_list = ['fitDiagnostics_partblind_uncorbtaglfhf_run2.root']
    for f in f_list:
        postfit = PostFit(f)
        save_pdf(f.replace('root','pdf'))(postfit.makeplots)() # decorator


class PostFit:
    '''
    script that makes post-fit and pre-fit 
    stack plots of data vs MC from fitDiagnostic file
    '''
    #
    #hists  = {}  # dictionary with format Xfit:ZhptX:process:values
    #edges  = {}  # dictionary with format ZhptX:edges
    fit_dir = 'fitdiag_roots/'
    gof_dir = 'fitdiag_roots/gof'

    def __init__(self,fitroo, kinem=None, t_labels=None):
        self.fitroo     = fitroo
        self.name       = '_'.join(fitroo.replace('.root','').split('_')[1:])
        self.gofroo     = f"{self.gof_dir}/higgsCombine_{self.name}.GoodnessOfFit.mH120.root"
        self.gofroo_toy = f"{self.gof_dir}/higgsCombine_{self.name}_toy.GoodnessOfFit.mH120.root"
        self.kinem      = kinem
        self.t_labels   = np.array(t_labels) if t_labels is not None else t_labels
        self.year       = re.search(r'(201\d|run2)',fitroo).group()
        self.lumi       = round(cfg.Lumi[self.year],1) if self.year != 'run2' else 137
        self.chi2_arr   = {'prefit':[], 'postfit':[]}
        self.hists = {}
        self.edges = {}
        self.load_info()

    def load_info(self):
        with uproot.open(self.fit_dir+self.fitroo) as roo:
            # store, prefit # shapes_prefit, shapes_fit_s
            def to_dict(p): # takes prefit or postfit
                pp_dict = {'prefit' :'shapes_prefit',
                           'postfit':'shapes_fit_s'}
                suf = (lambda b: b.decode().split(';')[0])
                self.hists[p] = {}
                for ch in roo[pp_dict[p]]:
                    ch_str = suf(ch)
                    self.hists[p][ch_str] = {}
                    kill_empty = (lambda h_: h_[roo[pp_dict[p]][ch]['total'].values > 0]) # default non nnqc post fit 
                    if ch_str not in self.edges and self.t_labels is None:
                        self.edges[ch_str]=np.array(roo[pp_dict[p]][ch]['total'].edges) # dummy use of TTBar
                        self.edges[ch_str]=np.append(self.edges[ch_str][0],self.edges[ch_str][1:][roo[pp_dict[p]][ch]['total'].values > 0]) # only kill empty right-most bin
                    elif self.t_labels is not None:
                        # address np.inf limit
                        t_labels = np.where(self.t_labels != np.inf, self.t_labels,
                                            2*self.t_labels[-2]-self.t_labels[-3])
                        self.edges[ch_str]=t_labels
                        kill_empty = (lambda h_: h_) # dont drop nnqc post fit bins
                    for hist in roo[pp_dict[p]][ch]:
                        hist_str = suf(hist)
                        h = roo[pp_dict[p]][ch][hist]
                        try:
                            self.hists[p][ch_str][hist_str] = {'values':kill_empty(h.values),  'err':kill_empty(np.sqrt(h.variances))}
                        except:
                            self.hists[p][ch_str][hist_str] = {'values':kill_empty(h.yvalues), 'errup':kill_empty(h.yerrorshigh),'errdw':kill_empty(h.yerrorslow)}
            #
            to_dict('prefit')
            to_dict('postfit')
        #
    #
    def makeplots(self, doPull=True, do1Pull=True,):
        self.doPull = doPull
        #self.do1Pull = do1Pull
        self.pulls = np.array([])
        for ch in self.hists['prefit']: # 'prefit is just a placeholder to loop over pt'
            self.init_axes(channel=ch)
            self.fig.suptitle(ch, fontsize=10)
            self.make_stackhist(ch)
            if doPull:
                self.init_axes(opt='pulls',channel=ch)
                self.fig.suptitle(ch, fontsize=10)
                self.make_pulls(ch)
        if self.doPull:
            self.plot_1dpulls()
            self.plot_1dgof()
        #print(f"data test statistic: prefit:{np.nansum(self.chi2_arr['prefit'])} | postfit: {np.nansum(self.chi2_arr['postfit'])}") 
        print(f"data test statistic: prefit:{-2*np.log(np.nanprod(self.chi2_arr['prefit']))} | postfit: {-2*np.log(np.nanprod(self.chi2_arr['postfit']))}") 

    def plot_1dpulls(self):
        fig, ax = plt.subplots()
        #print(self.pulls)
        self.pulls = self.pulls[(~np.isnan(self.pulls)) & (self.pulls != 0)]
        mu, std = norm.fit(self.pulls.flatten())
        #h = ax.hist(self.pulls.flatten(), 
        #            histtype='step',
        #            #bins=[-4,-3.5,-3.0,-2.5,-2.0,-1.25,-0.75,-0.25,
        #            #0.25,0.75,1.25,2.0,2.5,3,3.5,4])#[-4,-3,-2,-1,1,2,3,4])
        #            bins = 40,
        #            density=True)
        #xmin, xmax = ax.get_xlim()
        xmin, xmax = np.min(self.pulls), np.max(self.pulls)
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        ax.plot(x, p, 'k', linewidth=2)
        fig.text(x=.75,
                 y=.75,
                 s=f'Mean: {np.nanmean(self.pulls.flatten()):.3f}\nStd   : {np.nanstd(self.pulls.flatten()):.3f}')
        self.end_pullsgof(fig, ax, title='Pulls')
        #ax.set_xlabel('Pull'+('' if self.kinem is None else f' ({self.kinem})'))
        #ax.set_ylabel('Norm. / bin')
        #ax.xaxis.set_minor_locator(AutoMinorLocator())
        #ax.yaxis.set_minor_locator(AutoMinorLocator())
        #ax.grid(False)
        #fig.text(0.15,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 12)
        #fig.text(0.70,0.89, f'{self.lumi}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 12)
        #ax.set_title(f'Pulls in {self.year}')

        #fig.suptitle(f'Pulls in {self.year}')
        #plt.show()

    def plot_1dgof(self):
        fig, ax = plt.subplots()
        t     = uproot.open(self.gofroo)['limit']
        t_toy = uproot.open(self.gofroo_toy)['limit']
        gof_data =  t.array('limit')[0]
        gof_toy  =  t_toy.array('limit')
        mu, std = norm.fit(gof_toy.flatten())
        #xmin, xmax = ax.get_xlim()
        xmin, xmax = np.min(gof_toy), np.max(gof_toy)
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        ax.plot(x, p, 'k', linewidth=2)
        ax.hist(gof_toy.flatten(), 
                histtype='step',density=True)
        #print(gof_data)
        data_stat = gof_data# - mu)/std
        pscore    = len(gof_toy[gof_toy>=gof_data])/len(gof_toy)
        fig.text(
            x=.65,
            y=.75,
            s=f'Data Stat: {data_stat:.3f}\nP-score: {pscore:.3f}')
        ax.axvline(x=data_stat, color='red', linewidth='1', linestyle='-')
        self.end_pullsgof(fig, ax, title='GOF')

    def end_pullsgof(self, fig, ax, title='Pulls'):
        ax.set_xlabel(title+('' if self.kinem is None else f' ({self.kinem})'), fontsize=8)
        ax.set_ylabel('Norm. / bin')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(False)
        fig.text(0.15,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 10)
        fig.text(0.70,0.89, f'{self.lumi}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
        #ax.set_title(f'{title} for {self.year}')
    #
    def init_axes(self, opt='', channel=''):
        # init fig and four axes
        #fig, (ax1t, ax1b, ax2t, ax2b) = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig_ax_map = {
                '' : (lambda  : plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})), # 
            'pulls': (lambda  : plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})), # just 1 axes for pulls
            'blind' : plt.subplots
        }
        #fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig, axs = fig_ax_map[opt]()
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
        if opt == 'pulls':
            self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 10)
            self.fig.text(0.70,0.89, f'{lumi:.1f}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
        elif opt == 'blind':
            self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Simulation$", fontsize = 10)
            self.fig.text(0.70,0.89, f'{lumi:.1f}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 10)
        else:
            self.fig.text(0.085,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 8)
            self.fig.text(0.55,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 8)
            self.fig.text(0.30,0.89, f'{lumi:.1f}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 8)
            self.fig.text(0.765,0.89, f'{lumi:.1f}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 8)
        self.axs = np.array(axs) # [1:[t,b],2:[t,b]]
    #
    def do_stack(self, d, top_axs, ch):
        ycumm = None
        # stack portion
        #ordered_list = re.findall(rf'tt[H,Z]\d', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_2b','tt_bb','TTBar']
        ordered_list = re.findall(rf'tt[H,Z]\d', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_B','TTBar']
        #colors =  plt.cm.gist_rainbow(np.linspace(0,1,len(ordered_list)))
        colors =  plt.cm.tab20c(np.linspace(0,1,20))[:8]
        colors = np.append(colors, plt.cm.gist_rainbow(np.linspace(0,1,6)), axis=0)
        #
        #for j,k in enumerate(d):
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
            top_axs.fill_between(self.edges[ch],ycumm,ycumm-y,step="post", 
                                 linewidth=0, color=c, label=label)
            # add total error and data points
        top_axs.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=d['data']['values'],
                     xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,yerr=[d['data']['errdw'],d['data']['errup']], 
                     fmt='.', label='data', color='k')
        self.make_error_boxes(top_axs, (self.edges[ch][1:]+self.edges[ch][:-1])/2, d['total']['values'],
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err'], label='stat+sys')
        
    #
    def make_stackhist(self,ch='Zhpt2'):
        for i,p in zip(range(len(self.axs)),self.hists):
            axt = self.axs[0,i] # # 0,0 lt # 0,1 rt
            axb = self.axs[1,i]
            d = self.hists[p][ch]
            ## stack plot here
            self.do_stack(d, axt, ch)
            # --- for testing 
            #self.chi2_arr[p] = np.append(self.chi2_arr[p], np.power(d['data']['values']-d['total']['values'],2)/d['data']['values'])
            #self.chi2_arr[p] = np.append(self.chi2_arr[p], (np.power((d['data']['values']-d['total']['values']),2))/d['data']['values'])
            self.chi2_arr[p] = np.append(self.chi2_arr[p], np.exp( -1*np.power(d['data']['values']-d['total']['values'],2) / (2*(d['total']['values']-np.power(d['total']['err'],2))) ))
            # --- end testing 
            # add data/mc plot underneath
            y       =  d['data']['values']/d['total']['values']
            yerrdw  =  d['data']['errdw'] /d['total']['values']
            yerrup  =  d['data']['errup'] /d['total']['values']
            axb.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=y,
                         xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,yerr=[yerrdw,yerrup], 
                         fmt='.', label='data', color='k')
            self.make_error_boxes(axb, (self.edges[ch][1:]+self.edges[ch][:-1])/2, np.ones_like(d['total']['values']),
                                  xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err']/d['total']['values'], label='stat+sys')
            #
            self.endaxs(axt,axb,p, ch=ch)
            self.addlegend()
        #
        #plt.show()
    #
    def make_pulls(self,ch):
        # only for post fit
        axt = self.axs[0]
        axb = self.axs[1]
        d = self.hists['postfit'][ch]
        ## stack portion
        self.do_stack(d, axt ,ch)
        # Pull def = (Data - Pred.) / sqrt( Pred. + Pred.err^2 )
        #y =        (d['data']['values'] - d['total']['values']) / np.sqrt( d['total']['values'] + np.power( d['total']['err'] ,2))
        #print(d['data']['values'], d['total']['values'], np.power(d['total']['err'],2) )
        
        #print()
        #print(d['total']['values'])
        #print(d['total']['err'])
        #print()
        
        num = np.subtract(d['data']['values'], d['total']['values'])
        #den = np.sqrt( np.subtract( d['total']['values'], np.power( d['total']['err'] ,2) ))
        den = np.sqrt( abs(np.subtract( d['total']['values'], np.power( d['total']['err'] ,2) )))
        #y = (d['data']['values'] - d['total']['values']) / np.sqrt( d['total']['values'] - np.power( d['total']['err'] ,2))
        y = num/den
        #print(y)
        self.pulls = np.append(self.pulls, y)
        #yerrdw  =  (d['data']['errdw'] ) / np.sqrt( d['total']['values'] + np.power( d['total']['err'] ,2))
        #yerrup  =  (d['data']['errup'] ) / np.sqrt( d['total']['values'] + np.power( d['total']['err'] ,2))
        axb.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=y,
                     xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,#yerr=[yerrdw,yerrup], 
                     fmt='.', label='data', color='k')
        #self.make_error_boxes(axb, (self.edges[ch][1:]+self.edges[ch][:-1])/2, np.ones_like(d['total']['values']),
        #                      xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err']/d['total']['values'], label='stat+sys')
        self.endaxs(axt,axb,'postfit',self.doPull, ch)
        self.addlegend(opt='pulls')
        

    def endaxs(self,axt,axb,p,doPull=False, ch='Zhpt1'):
        #axt
        axt.set_xlim(0,self.edges[ch][-1])
        axt.set_ylim(0.)

        #if axt.get_ylim()[1] > 1000:
        axt.set_yscale('log')
        axt.set_ylim(ymin=1., ymax=axt.get_ylim()[1] * 10)
        #else:
        #     axt.yaxis.set_minor_locator(AutoMinorLocator())

        axt.set_ylabel('Events / bin', fontsize=8)
        axt.set_title(rf'${p}$', y=1.0, pad=-14, fontsize=8)
        axt.xaxis.set_minor_locator(AutoMinorLocator(5))

        axt.tick_params(which='both', direction='in', top=True, right=True)
        #axb

        axb.xaxis.set_minor_locator(AutoMinorLocator(5))
        axb.yaxis.set_minor_locator(AutoMinorLocator())
        axb.tick_params(which='both', direction='in', top=True, right=True)
        if not doPull:
            axb.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
            axb.set_ylim(0.5,1.5)
            axb.yaxis.set_major_locator(FixedLocator([.75,1,1.25]))
            axb.yaxis.set_major_formatter(FormatStrFormatter('%g'))
            axb.set_ylabel('data/MC' if 'post' not in p else 'data/pred.', fontsize=8)
        else:
            axb.axhline(0, color='red', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
            axb.set_ylim(-3.0,3.0)
            axb.yaxis.set_major_locator(FixedLocator([-2,-1,0.0,1,2]))
            axb.yaxis.set_major_formatter(FormatStrFormatter('%g'))
            axb.set_ylabel('Pull', fontsize=8)
        axb.set_xlabel(('' if self.kinem is None else self.kinem+' ')+('bin #' if self.t_labels is None else ''), fontsize=8)
        axb.grid(True)
        
    def addlegend(self, opt=''):
        leg_map = {
            ''     : (lambda : self.axs[0,1]), # defualt
            'pulls': (lambda : self.axs[0]), # pulls
            'blind': (lambda : self.axs), # pulls
            }
        #ax = self.axs[0,1]
        ax = leg_map[opt]()
        handles, labels = ax.get_legend_handles_labels()
        hatch_patch = Patch(hatch=10*'X', label='stat+sys',  fc='w')
        handles = handles + [hatch_patch]
        labels  = labels + ['stat+sys.']
        ax.legend(handles,labels, bbox_to_anchor=(1.00,1), 
                  fontsize='xx-small', framealpha = 0, loc='upper left')


    @staticmethod
    def make_error_boxes(ax, xdata, ydata, xerror, yerror,  facecolor='r',
                         edgecolor='None', alpha=0.0, hatch=10*'X', label=''):
        errorboxes = []
        for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
            rect = Rectangle((x - xe, y - ye), 2*xe, 2*ye)
            errorboxes.append(rect)
        pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                             edgecolor=edgecolor, hatch=hatch, label=label, zorder=1.5)
        # Add collection to axes
        ax.add_collection(pc)

if __name__ == '__main__':
    #
    main()
