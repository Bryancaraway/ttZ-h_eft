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

dcdir = f'{sys.path[1]}/Higgs-Combine-Tool/'
sys_to_plot = ['jesRelativeSample{y}','jesHF{y}','jesAbsolute{y}','jesEC2{y}',
               'jesBBEC1{y}', 'jms_{y}', 'jmr_{y}','ak4jer_{y}', 'ak8jer_{y}',
               'jesRelativeBal', 'jesFlavorQCD', 'jesHF', 'jesAbsolute', 'jesEC2', 
               'jesBBEC1', 'fsr', 'UE', 'hdamp']
sample = 'TTBar'
#sample = 'ttZ3'
year   = '2017'

@save_pdf(f'qcsmooth_before_after_{sample}_{year}.pdf')
def main():
    # can set qcsmooth to 'qc' or 'smooth'
    qcsmooth = PlotQCSmooth(sample,'Zhpt1',year, qcsmooth='smooth')
    qcsmooth.makeplots()
    qcsmooth = PlotQCSmooth(sample,'Zhpt2',year, qcsmooth='smooth')
    qcsmooth.makeplots()
    qcsmooth = PlotQCSmooth(sample,'Zhpt3',year, qcsmooth='smooth')
    qcsmooth.makeplots()

    

class PlotQCSmooth:

    def __init__(self, process, channel, year, qcsmooth='smooth'):
        self.process = process
        self.channel = channel
        self.year    = year
        self.roofiles  = [f'{dcdir}/datacard_{opt}{qcsmooth}_{year}.root' for opt in ['no','']]
        self.hists     = {}
        self.edges     =  None
        self.sys_to_plot = [sys.replace('{y}',f'{self.year}')  for sys in sys_to_plot]
        self.sys_names = []
        self.hist_dict = {}
        self.load_hists()

    def load_hists(self):
        for roofile in self.roofiles:
            self.hists_getter(roofile)
            self.hist_dict['before' if 'no' in roofile else 'after'] = self.hists 
            self.hists = {} # clear hist dict

    def hists_getter(self,roofile):
        with uproot.open(roofile) as roo:
            pre = f'{self.channel}_{self.process}'
            sys_names = [k.decode().replace('Up','').replace('Down','') for k in roo.keys()]
            self.sys_names = list(dict.fromkeys(re.findall(rf'{pre}_\w*' ,' '.join(sys_names))))
            self.hists['nom'] = np.array(roo[pre].values)   # sumw
            self.hists['err'] = (np.sqrt(roo[pre].variances)+self.hists['nom'])/self.hists['nom'] # sumw2
            if self.edges is None : self.edges = np.array(roo[pre].edges)
            #for sys in self.sys_to_plot:
            for sys in self.sys_names:
                for ud in ['Up','Down']:
                    #try:
                    #self.hists[sys+ud] = np.array(roo[pre+f'_{sys}'+ud].values)/self.hists['nom']
                    self.hists[sys+ud] = np.array(roo[sys+ud].values)/self.hists['nom']
                    if 'hdamp' in sys or 'UE' in sys:
                        #self.hists[sys+ud+'_err'] = (np.sqrt(roo[pre+f'_{sys}'+ud].variances)+np.array(roo[pre+f'_{sys}'+ud].values))/np.array(roo[pre+f'_{sys}'+ud].values)
                        self.hists[sys+ud+'_err'] = (np.sqrt(roo[sys+ud].variances)+np.array(roo[sys+ud].values))/np.array(roo[sys+ud].values)
                    #except KeyError: pass # for sys not in root file

    def makeplots(self):
        #for sys in self.sys_to_plot:
        for sys in self.sys_names:
            try:
                fig, axs = self.initPlt()
                for ax,sm in zip(axs,self.hist_dict):
                    self.make_step(sys,ax,sm)
                self.endPlt(fig,axs)
            except KeyError: pass # for sys not in root file    

    def make_step(self, sys, ax, opt):
        ax.set_title(rf'${opt}$', y=1.0, pad=-14, fontsize=10)
        ax.set_xlabel('Bin number')
        ax.axhline(1, color='r', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        if 'hdamp' in sys or 'UE' in sys:
            for ud,c in zip(['Up','Down'],['tab:blue','tab:orange']):
                self.make_error_boxes(ax, (self.edges[1:]+self.edges[:-1])/2, self.hist_dict[opt][sys+ud],
                                      (self.edges[1:]-self.edges[:-1])/2, 
                                      abs(self.hist_dict[opt][sys+ud+'_err']-self.hist_dict[opt][sys+ud]), 
                                      label='stat unc.', facecolor=c)
        else:
            self.make_error_boxes(ax, (self.edges[1:]+self.edges[:-1])/2, self.hist_dict[opt]['nom']/self.hist_dict[opt]['nom'],
                                  (self.edges[1:]-self.edges[:-1])/2, abs(self.hist_dict[opt]['err']-1), label='stat unc.')
        
        for ud in ['Up','Down']:
            ax.step(
                x = self.edges,
                y = np.append(self.hist_dict[opt][sys+ud],0),
                where='post', label=f'{sys}{ud}'
            )
        ax.set_xlim(self.edges[0],self.edges[-1])
        ax.set_ylim(0.0,2.0)
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        handles, labels = ax.get_legend_handles_labels()
        hatch_patch = Patch(hatch=10*'X', label='stat unc.',  fc='w')
        handles = handles + [hatch_patch]
        labels  = labels + ['stat unc.']
        ax.legend(handles,labels, #bbox_to_anchor=(1.00,1), 
                  fontsize=8, loc='lower left')
        #ax.legend(loc=8,fontsize=10)
        ax.grid(True)

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

            
    def initPlt(self):
        fig, axs = plt.subplots(1,2, sharey=True)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.3,
            wspace=0.02
        )
        return fig, axs
        
    def endPlt(self, fig, axs):
        fig.suptitle(f'Smoothing for {self.process}, {self.channel}({self.year})')
        axs[0].set_ylabel('Variation/Nom')




if __name__ == '__main__':
    main()

