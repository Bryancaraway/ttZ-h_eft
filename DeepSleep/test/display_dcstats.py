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
dc_shapes = 'datacard_year.root'
ch_name = 'Zhpt3'
#channels = ['Zhpt2','Zhpt3']
channels = ['Zhpt3']

ttZ_list = [f'ttZ{i}' for i in range(4)]
ttH_list = [f'ttH{i}' for i in range(4)]
#sample_list = ['tt_B','TTBar','ttX','single_t','VJets'] + ttZ_list + ttH_list
#sample_list = ['tt_B','TTBar','ttX','single_t'] + ttZ_list + ttH_list
sample_list = ttZ_list  + ttH_list + ['single_t']


#@save_pdf(f'all_mcstats.pdf')
@save_pdf(f'singlet_mcstats.pdf')
def main():
    for ch in channels:
        for y in cfg.Years:
            mcstat = Plot_mcstats(y, ch)
            mcstat.makeplots()
            for bin_i in [21,22]:
                pie_chart = Plot_pie_mcstats(y, ch, bin_i)
                #pie_chart.makeplots()               

class Plot_mcstats:
    
    def __init__(self, year, ch):
        self.roofile = dc_dir + dc_shapes.replace('year',year)
        self.year    = year
        self.channel = ch
        self.hist_errs  = {} # initiate 
        self.hist_vals  = {}
        self.hist_edges = []
        self.load_hists()
    
    def load_hists(self):
        with uproot.open(self.roofile) as roo:
            for sample in sample_list: 
                hist_name = f'{self.channel}_{sample}'
                self.hist_vals[sample] = roo[hist_name].values
                self.hist_errs[sample] = roo[hist_name].variances
            self.hist_errs['ttZ'] = sum([self.hist_errs[ttZ] for ttZ in  ttZ_list if ttZ in sample_list])
            self.hist_vals['ttZ'] = sum([self.hist_vals[ttZ] for ttZ in  ttZ_list if ttZ in sample_list])
            self.hist_errs['ttH'] = sum([self.hist_errs[ttH] for ttH in  ttH_list if ttH in sample_list])
            self.hist_vals['ttH'] = sum([self.hist_vals[ttH] for ttH in  ttH_list if ttH in sample_list])
            for sig in ttZ_list+ttH_list:
                if sig in self.hist_errs:
                    del self.hist_errs[sig], self.hist_vals[sig]
            self.hist_edges = roo[hist_name].edges # hist_name in local scope?


    def makeplots(self):
        fig, ax = self.initPlt()
        #
        ch_edges = np.arange(66+1-len(self.hist_edges),66+1)
        color_dict = {'ttZ':'tab:blue',
                      'ttH':'tab:orange',
                      'single_t':'tab:brown'}
        for sample in self.hist_errs:
            #if sample == 'ttZ':
            #    print(sample, self.year)
            #    for_olaf = np.array([f'{i:.3f} +/- {j:.3f}' for i,j in zip(self.hist_vals[sample],np.sqrt(self.hist_errs[sample]))])[[13,17,21]]
            #    print(f'bin 14 {for_olaf[0]} | bin 18 {for_olaf[1]} | bin 22: {for_olaf[2]}')
            if sample == 'ttH':
                print(sample, self.year)
                #for_olaf = np.array([f'{i:.3f} +/- {j:.3f}' for i,j in zip(self.hist_vals[sample],np.sqrt(self.hist_errs[sample]))])[[14,18,22]]
                #for_olaf = np.round(np.sum(np.array(self.hist_vals['ttZ'])[[13,17,21]])/np.sum(np.array(self.hist_vals[sample])[[14,18,22]]),3)
                #for_olaf = np.round(np.sum(np.array(self.hist_vals['ttZ'])[12:])/np.sum(np.array(self.hist_vals[sample])[12:]),3)
                #for_olaf_err = np.round(
                #    for_olaf* np.sqrt(
                #        np.power(np.array(self.hist_errs['ttZ'])[[13,17,21]]/np.array(self.hist_vals['ttZ'])[[13,17,21]],2) + \
                #        np.power(np.array(self.hist_errs['ttH'])[[14,18,22]]/np.array(self.hist_vals['ttH'])[[14,18,22]],2) 
                #    ), 3)
                #print(for_olaf)
                #print(f'bin 15 {for_olaf[0]} | bin 19 {for_olaf[1]} | bin 23: {for_olaf[2]}')
                #print(f'bin 15 {for_olaf[0]} +/- {for_olaf_err[0]} | bin 19 {for_olaf[1]} +/- {for_olaf_err[1]} | bin 23: {for_olaf[2]} +/- {for_olaf_err[2]}')
                
            ax.step(
                x = ch_edges,
                #y = np.append(np.sqrt(self.hist_errs[sample]),0),
                #y = np.append(100*np.sqrt(self.hist_errs[sample])/self.hist_vals[sample],0),
                y = np.append(self.hist_vals[sample],0),
                color = color_dict.get(sample,'k'),
                where='post', label=sample,
            )
            ax.errorbar(
                x = (ch_edges[1:]+ch_edges[:-1])/2,

                #y = np.append(100*np.sqrt(self.hist_errs[sample])/self.hist_vals[sample],0),
                y = self.hist_vals[sample],
                yerr = np.sqrt(self.hist_errs[sample]),
                fmt='none', 
                color = color_dict.get(sample,'k'),
                #where='post', label=sample,
            )
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(ch_edges[0],ch_edges[-1])
        #ax.set_ylim(0,105)
        #ax.set_yscale('log')
        ax.legend()
        ax.grid(True)

        #
        self.endPlt(fig,ax)
        #plt.show()
        
        

        
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
        fig.suptitle(f'MC stat in channel {self.channel}, {self.year}')
        ax.set_ylabel('MC yield / bin')
        #ax.set_ylabel(r'$\frac{\mathrm{|MC\;stat\;unc.|}}{\mathrm{MC\;yield}} \times 100$')

class Plot_pie_mcstats(Plot_mcstats):
    def __init__(self, year, ch, bin_n):
        self.bin_n = bin_n
        super().__init__(year,ch)

    def makeplots(self):
        fig, ax = self.initPlt()
        ax.pie(np.sqrt([self.hist_errs[sample][self.bin_n] for sample in self.hist_errs]), labels=self.hist_errs.keys(), normalize=True, autopct='%1.1f%%')
        ax.axis('equal')
        fig.suptitle(f'MC stat error in channel : {self.channel}, bin : {self.bin_n}, year :  {self.year}')
        #plt.show()

if __name__ == '__main__':
    main()
