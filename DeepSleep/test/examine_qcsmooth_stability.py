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
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            

sys_to_plot = ['hdamp','UE']

#sample = 'TTBar'
#year = '2016'

# smoothing means sum(up/down) / sum(nom)
@save_pdf('smoothing_stability.pdf')
def main():
    for sample in ['TTBar','tt_B']:
        for year in cfg.Years:
            for ch in ['Zhpt1','Zhpt2','Zhpt3']:
                smooth_stable = SmoothStability(sample,ch,year,qcsmooth='smooth')
                smooth_stable.makeplots()

class SmoothStability(PlotQCSmooth):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hist_dict = self.hist_dict['before'] # dont care about after

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
            self.hists['nom'] = np.array(roo[pre].values)   # sumw
            self.hists['err'] = np.array(roo[pre].variances)
            if self.edges is None : self.edges = np.array(roo[pre].edges)
            #for sys in self.sys_to_plot:
            for sys in self.sys_names:
                for ud in ['Up','Down']:
                    #try:
                    self.hists[sys+ud] = np.array(roo[sys+ud].values)
                    if 'hdamp' in sys or 'UE' in sys:
                        self.hists[sys+ud+'_err'] = np.array(roo[sys+ud].variances)
                    #except KeyError: pass # for sys not in root file

    def makeplots(self):
        for sys in self.sys_names:
            self.throw_toys(sys)

    def throw_toys(self, sys):
        nom = np.nansum(self.hist_dict['nom'])
        fig, ax = plt.subplots()
        for ud_str in ['Up','Down']:
            ud = self.hist_dict[sys+ud_str]
            err = np.sqrt(self.hist_dict[sys+ud_str+'_err'])
            #print(ud)
            #print(err)        
            toys = np.random.normal(ud,err, size = (10000,len(ud)))
            smooth_toys = np.nansum(toys, axis=1)
            #print(smooth_toys, smooth_toys.shape)
            #print(np.mean(smooth_toys), np.std(smooth_toys), nom, np.nansum(ud))
            ax.hist(smooth_toys/nom, range=(0.90,1.10), bins=100, histtype='step', label=f'{sys}{ud_str} (10K toys)')

        ax.axvline(1, color='red', linewidth='1', linestyle='-', snap=True, alpha=0.5)
        ax.set_ylabel('Count / bin')
        ax.set_xlabel('Toy Variation / Nominal')
        fig.suptitle(f'{self.channel}_{self.year}')
        ax.legend(fontsize=8)
        #plt.show()
        #exit()
        



if __name__ == '__main__':
    main()
