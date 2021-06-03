import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')

import numpy as np
np.random.seed(1)
import re
import uproot
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

dc_dir = f'{sys.path[1]}/Higgs-Combine-Tool/'
dc_shapes = 'datacard_2018.root'
ch_name = 'Zhpt3'

ttZ_list = [f'ttZ{i}' for i in range(4)]
ttH_list = [f'ttH{i}' for i in range(4)]
sample_list = ['tt_B','TTBar','ttX','single_t','VJets'] + ttZ_list + ttH_list



def main():
    mcstat = Plot_mcstats()
    save_pdf(f'{ch_name}_mcstats.pdf')(mcstat.makeplots)()

class Plot_mcstats:
    
    def __init__(self):
        self.roofile = dc_dir + dc_shapes
        self.hist_errs  = {} # initiate 
        self.hist_edges = []
        self.load_hists()
    
    def load_hists(self):
        with uproot.open(self.roofile) as roo:
            for sample in sample_list: 
                hist_name = f'{ch_name}_{sample}'
                self.hist_errs[sample] = roo[hist_name].variances
            self.hist_errs['ttZ'] = sum([self.hist_errs[ttZ] for ttZ in  ttZ_list])
            self.hist_errs['ttH'] = sum([self.hist_errs[ttH] for ttH in  ttH_list])
            for sig in ttZ_list+ttH_list:
                del self.hist_errs[sig]
            self.hist_edges = roo[hist_name].edges # hist_name in local scope?


    def makeplots(self):
        fig, ax = self.initPlt()
        #
        for sample in self.hist_errs:
            ax.step(
                x = self.hist_edges,
                y = np.append(np.sqrt(self.hist_errs[sample]),0),
                where='post', label=sample,
            )
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(self.hist_edges[0],self.hist_edges[-1])
        #ax.set_yscale('log')
        ax.legend()
        ax.grid(True)
        #
        self.endPlt(fig,ax)
        

        
    def initPlt(self):
        fig, ax = plt.subplots()
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.3,
            wspace=0.02
        )
        return fig, ax
        
    def endPlt(self, fig, ax):
        fig.suptitle(f'MC stat in channel {ch_name}')
        ax.set_ylabel('|MC stat unc.|')

if __name__ == '__main__':
    main()
