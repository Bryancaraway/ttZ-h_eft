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
from lib.fun_library import save_pdf
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            


def main():
    # fitDiagnostics_test_2016.root, fitDiagnostics_test_2017.root, fitDiagnostics_test_2018.root, sfitDiagnostics_test_run2.root
    #PostFit("fitDiagnostics_test_2016.root").makeplots()
    f_list = [f"fitDiagnostics_test_{i}.root" for i in ['2016','2017','2018','run2']]
    for f in f_list:
        postfit = PostFit(f)
        save_pdf(f.replace('root','pdf'))(postfit.makeplots)() # decorator


class PostFit:
    '''
    script that makes post-fit and pre-fit 
    stack plots of data vs MC from fitDiagnostic file
    '''
    #
    hists  = {}  # dictionary with format Xfit:ZhptX:process:values
    edges  = {}  # dictionary with format ZhptX:edges

    def __init__(self,fitroo):
        self.fitroo = fitroo
        self.year = re.search(r'(201\d|run2)',fitroo).group()
        self.lumi = round(cfg.Lumi[self.year],1) if self.year != 'run2' else 137
        self.load_info()

    def load_info(self):
        with uproot.open(self.fitroo) as roo:
            # store, prefit # shapes_prefit, shapes_fit_s
            def to_dict(p): # takes prefit or postfit
                pp_dict = {'prefit' :'shapes_prefit',
                           'postfit':'shapes_fit_s'}
                suf = (lambda b: b.decode().split(';')[0])
                self.hists[p] = {}
                for pt in roo[pp_dict[p]]:
                    pt_str = suf(pt)
                    self.hists[p][pt_str] = {}
                    if pt_str not in self.edges:
                        self.edges[pt_str]=roo[pp_dict[p]][pt]['TTBar'].edges
                    for hist in roo[pp_dict[p]][pt]:
                        hist_str = suf(hist)
                        h = roo[pp_dict[p]][pt][hist]
                        try:
                            self.hists[p][pt_str][hist_str] = {'values':h.values,  'err':np.sqrt(h.variances)}
                        except:
                            self.hists[p][pt_str][hist_str] = {'values':h.yvalues, 'errup':h.yerrorshigh,'errdw':h.yerrorslow}
            #
            to_dict('prefit')
            to_dict('postfit')
        #
    #
    def makeplots(self):
        for pt in self.hists['prefit']: # 'prefit is just a placeholder to loop over pt'
            self.init_axes()
            self.fig.suptitle(pt)
            self.make_stackhist(pt)
    #
    def init_axes(self):
        # init fig and four axes
        #fig, (ax1t, ax1b, ax2t, ax2b) = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.0,
            wspace=0.3
        )
        self.fig = fig
        self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Preliminary$", fontsize = 12)
        self.fig.text(0.635,0.89, f'{self.lumi}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 12)
        self.axs = np.array(axs) # [1:[t,b],2:[t,b]]
    #
    def make_stackhist(self,pt='Zhpt2'):
        for i,p in zip(range(len(self.axs)),self.hists):
            axt = self.axs[0,i] # # 0,0 lt # 0,1 rt
            axb = self.axs[1,i]
            d = self.hists[p][pt]
            colors =  plt.cm.gist_rainbow(np.linspace(0,1,len(d)))
            ycumm = None
            # stack portion
            for j,k in enumerate(d):
                if 'data' in k or 'total' in k: continue
                y = np.append(d[k]['values'],0)
                if ycumm is not None:
                    ycumm += y
                else:
                    ycumm = y 
                axt.fill_between(self.edges[pt],ycumm,ycumm-y,step="post", 
                                 linewidth=0, color=colors[j], label=k)
            # add total error and data points
            axt.errorbar(x=(self.edges[pt][1:]+self.edges[pt][:-1])/2, y=d['data']['values'],
                         xerr=(self.edges[pt][1:]-self.edges[pt][:-1])/2 ,yerr=[d['data']['errdw'],d['data']['errup']], 
                         fmt='.', label='data', color='k')
            self.make_error_boxes(axt, (self.edges[pt][1:]+self.edges[pt][:-1])/2, d['total']['values'],
                                  xerror=(self.edges[pt][1:]-self.edges[pt][:-1])/2, yerror=d['total']['err'], label='stat+sys')
            # add data/mc plot underneath
            y    = d['data']['values']/d['total']['values']
            #yerrdw = y*np.sqrt( np.power(d['data']['errdw']/d['data']['values'],2) + np.power(d['total']['err']/d['total']['values'],2) )
            #yerrup = y*np.sqrt( np.power(d['data']['errup']/d['data']['values'],2) + np.power(d['total']['err']/d['total']['values'],2) )
            yerrdw  =  d['data']['errdw']/d['total']['values']
            yerrup  =  d['data']['errup']/d['total']['values']
            axb.errorbar(x=(self.edges[pt][1:]+self.edges[pt][:-1])/2, y=y,
                         xerr=(self.edges[pt][1:]-self.edges[pt][:-1])/2 ,yerr=[yerrdw,yerrup], 
                         fmt='.', label='data', color='k')
            self.make_error_boxes(axb, (self.edges[pt][1:]+self.edges[pt][:-1])/2, np.ones_like(d['total']['values']),
                                  xerror=(self.edges[pt][1:]-self.edges[pt][:-1])/2, yerror=d['total']['err']/d['total']['values'], label='stat+sys')
            #
            self.endaxs(axt,axb,p)
            self.addlegend()
        #plt.show()

    def endaxs(self,axt,axb,p):
        #axt
        axt.set_xlim(0,16)
        axt.set_ylim(0)
        axt.set_ylabel('Events / bin')
        axt.set_title(rf'${p}$', y=1.0, pad=-14)
        axt.xaxis.set_minor_locator(AutoMinorLocator())
        axt.yaxis.set_minor_locator(AutoMinorLocator())
        axt.tick_params(which='both', direction='in', top=True, right=True)
        #axb
        axb.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        axb.xaxis.set_minor_locator(AutoMinorLocator())
        axb.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        axb.yaxis.set_major_locator(FixedLocator([.75,1,1.25]))
        axb.yaxis.set_minor_locator(AutoMinorLocator())
        axb.tick_params(which='both', direction='in', top=True, right=True)
        axb.set_ylim(0.5,1.5)
        axb.set_ylabel('data/MC' if 'post' not in p else 'data/pred.')
        axb.set_xlabel('bin #')
        axb.grid(True)
        
    def addlegend(self):
        ax = self.axs[0,1]
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
