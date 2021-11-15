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
from display_qcsmooth import PlotQCSmooth

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib import transforms
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            

sys_to_plot = ['hdamp','UE']
channels = ['Zhpt1','Zhpt2','Zhpt3']
#sample = 'TTBar
#year = '2016'

# smoothing means sum(up/down) / sum(nom)
#@save_pdf('smoothing_stability.pdf')
@save_pdf("dedsys_plots_TOP.pdf")
def main():
    for sample in ['TTBar','tt_B']:
        final_hists = None
        for year in cfg.Years:
            hists = {}
            hists_toys = {}
            for ch in ['Zhpt1','Zhpt2','Zhpt3']:
                smooth_stable = SmoothStability(sample,ch,year,qcsmooth='smooth')
                #smooth_stable.makeplots()
                hists.update(smooth_stable.hist_dict)
                hists_toys.update(smooth_stable.hist_toy_dict)
            merged_hists = make_plots_for_TOP(hists, hists_toys, sample,year) # pass to fuction which will make plots
            if final_hists is None:
                final_hists = merged_hists
            else:
                final_hists = {k: final_hists[k]+merged_hists[k] for k in final_hists}
            #print(final_hists['nom_all'])
        print(f"Run 2 rate uncertainty for {sample}")
        for sys in ['hdamp','UE','hdamp_ttbb']:
            if f'{sample}_{sys}Up_all' in final_hists:
                up = sum(final_hists[f'{sample}_{sys}Up_all'])/sum(final_hists['nom_all'])
                dn = sum(final_hists[f'{sample}_{sys}Down_all'])/sum(final_hists['nom_all'])
                print(f'{sys} Up/Down : \t', round(up/np.sqrt(up*dn),3),'/' , round(dn/np.sqrt(up*dn),3))
            #exit()
                

class SmoothStability(PlotQCSmooth):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hist_dict = self.hist_dict['before'] # dont care about after
        self.maketoys()

    def load_hists(self):
        for roofile in [self.roofiles[0]]: # only get no smooth plots
            self.hists_getter(roofile)
            self.hist_dict['before' if 'no' in roofile else 'after'] = self.hists 
            self.hists = {} # clear hist dict

    def hists_getter(self,roofile):
        with uproot.open(roofile) as roo:
            pre = f'{self.channel}_{self.process}'
            sys_names = [k.decode().replace('Up','').replace('Down','') for k in roo.keys()]
            self.sys_names = list(dict.fromkeys(re.findall(rf'{pre}_\w*' ,' '.join(sys_names))))
            self.sys_names = [name for name in self.sys_names for sys in sys_to_plot if sys in name]
            self.hists[f'{self.channel}_nom'] = np.array(roo[pre].values)   # sumw
            self.hists[f'{self.channel}_err'] = np.array(roo[pre].variances)
            if self.edges is None : self.edges = np.array(roo[pre].edges)
            #for sys in self.sys_to_plot:
            for sys in self.sys_names:
                for ud in ['Up','Down']:
                    #try:
                    self.hists[sys+ud] = np.array(roo[sys+ud].values)
                    if 'hdamp' in sys or 'UE' in sys:
                        self.hists[sys+ud+'_err'] = np.array(roo[sys+ud].variances)
                    #except KeyError: pass # for sys not in root file

    def maketoys(self):
        self.hist_toy_dict = {}
        for sys in self.sys_names:
            self.throw_toys(sys)

    def makeplots(self):
        for sys in self.sys_names:
            nom = np.nansum(self.hist_dict[f'{self.channel}_nom'])
            fig, ax = plt.subplots()
            for ud_str in ['Up','Down']:
                ax.hist(self.hist_toy_dict[f'{self.channel}_{sys}{ud_str}']/nom, range=(0.90,1.10), 
                        bins=100, histtype='step', label=f'{sys}{ud_str} (10K toys)')
            ax.axvline(1, color='red', linewidth='1', linestyle='-', snap=True, alpha=0.5)
            ax.set_ylabel('Count / bin')
            ax.set_xlabel('Toy Variation / Nominal')
            fig.suptitle(f'{self.channel}_{self.year}')
            ax.legend(fontsize=8)
    

    def throw_toys(self, sys):
        for ud_str in ['Up','Down']:
            ud = self.hist_dict[sys+ud_str]
            err = np.sqrt(self.hist_dict[sys+ud_str+'_err'])
            toys = np.random.normal(ud,err, size = (10000,len(ud)))
            smooth_toys = np.nansum(toys, axis=1)
            self.hist_toy_dict[f'{sys}{ud_str}'] = smooth_toys
            #self.hist_toy_dict[f'{self.channel}_{sys}{ud_str}_mean'] = np.mean(smooth_toys)
            #self.hist_toy_dict[f'{self.channel}_{sys}{ud_str}_std2'] = np.power(np.std(smooth_toys),2)
        

def make_plots_for_TOP(hists, hists_toys, sample, year):
    hists = combine_hists(hists)
    #hists_toys = combine_toy_hists(hists_toys)
    #for k in hists_toys:
    #    if 'err' in k:
    #        print(k)
    #        print(np.sqrt(hists[k]))
    #        print(hists_toys[k])
    #exit()
    pt_bins = np.array([200,300,450,600]) # 600 includes overflow
    rate_var = {
        'hdamp':[1.049,0.954], # Up,Down
        'hdamp_ttbb': [1.027,0.974],
        'UE': [1.01,0.99],
    }
    
    for sys in sys_to_plot:
        if sys == 'UE' and sample == 'tt_B': continue
        sys_name = sys if sample != 'tt_B' else sys+'_ttbb'
        # plot all unmerged
        fig, ax = plt.subplots()
        for i,ud_str in enumerate(['Up','Down']):
            x_edges = np.arange(len(hists['nom_all'])+1)
            ax.errorbar(x=.1*i+(x_edges[1:]+x_edges[:-1])/2, y=hists[f'{sample}_{sys_name}{ud_str}_all']/hists['nom_all'],
                        xerr=(x_edges[1:]-x_edges[:-1])/2,
                        yerr=np.sqrt(hists[f'{sample}_{sys_name}{ud_str}_err_all'])/hists['nom_all'],
                        fmt='.', label=f'{sys_name}{ud_str}_all')
        #
        trans = ax.transAxes + transforms.ScaledTranslation(0/72, 3/72, fig.dpi_scale_trans)
        ax.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        ax.text(0.06, 0.70, 'Z/H pT 200-300', transform=trans)
        ax.axvline(18, color='k', linewidth='1', linestyle='-', snap=True)
        ax.text(0.37, 0.70, 'Z/H pT 300-450', transform=trans)
        ax.axvline(42, color='k', linewidth='1', linestyle='-',  snap=True)
        ax.text(0.74, 0.70, 'Z/H pT 450-inf', transform=trans)
        ax.set_xlabel('template bins')
        ax.set_ylabel('Variation / Nominal')
        ax.set_xlim(-1,67)
        ax.set_ylim(0,5)
        fig.suptitle(f'{sample}_{year}')
        ax.grid(True)
        plt.legend()
        #plt.show()
        #exit()
        # plot merged
        fig, ax = plt.subplots()
        for i,ud_str in enumerate(['Up','Down']):
            ax.errorbar(x=3*i+(pt_bins[1:]+pt_bins[:-1])/2, y=hists[f'{sample}_{sys_name}{ud_str}']/hists['nom'],
                        xerr=(pt_bins[1:]-pt_bins[:-1])/2,
                        yerr=np.sqrt(hists[f'{sample}_{sys_name}{ud_str}_err'])/hists['nom'],
                        zorder=3,
                        fmt='.', label=f'{sys_name}{ud_str}')
            #ax.errorbar(x=3*i+(pt_bins[0]+pt_bins[-1])/2, y=sum(hists[f'{sample}_{sys_name}{ud_str}'])/sum(hists['nom']),
            #            xerr=(pt_bins[-1]-pt_bins[0])/2,
            #            yerr=np.sqrt(sum(hists[f'{sample}_{sys_name}{ud_str}_err']))/sum(hists['nom']),
            #            fmt='s', alpha=0.75, label=f'{sys_name}{ud_str}_rate')
        print(sys)
        ax.axhline(rate_var[sys][0], color='red', linewidth='1.5', linestyle='-', snap=True, label=f'{sys_name}_rate')
        ax.axhline(rate_var[sys][1], color='red', linewidth='1.5', linestyle='-', snap=True, )#label=f'{sys_name}Down_rate')
        ax.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        ax.set_xlabel('Z/H pT (GeV)')
        ax.set_ylabel('Variation / Nominal')
        ax.set_xlim(200,600)
        ax.set_ylim(min(0.85,ax.get_ylim()[0]),max(1.15,ax.get_ylim()[1]))
        #ax.set_yscale('log')
        fig.suptitle(f'{sample}_{year}')
        ax.grid(True)
        plt.legend()
        #plt.show()
        #exit()
    return hists
    


def combine_hists(hists):
    base_keys = [key.replace(rf'{channels[0]}_', '') for key in re.findall(rf'{channels[0]}_\w*', ' '.join(hists.keys()))]
    out_hist = {}
    for bkey in base_keys:
        out_hist[bkey] = np.array([np.nansum(hists[f'{ch}_'+bkey]) for ch in channels])
        out_hist[bkey+'_all'] = np.concatenate([hists[f'{ch}_'+bkey] for ch in channels])
    return out_hist
    
def combine_toy_hists(hists):
    base_keys = [key.replace(rf'{channels[0]}_', '') for key in re.findall(rf'{channels[0]}_\w*', ' '.join(hists.keys()))]
    out_hist = {}
    for bkey in base_keys:
        out_hist[bkey] = np.array([np.mean(hists[f'{ch}_'+bkey]) for ch in channels])
        out_hist[bkey+'_err'] = np.array([np.std(hists[f'{ch}_'+bkey]) for ch in channels])
    return out_hist

if __name__ == '__main__':
    main()
