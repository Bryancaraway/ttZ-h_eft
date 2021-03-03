import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import re
from lib.fun_library import save_pdf
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg

'''
Script that demonstrates how wc impacts
our expected signal in terms of pt distribution
'''


class main():
    pt_bins = {
        0 : [0,200],
        1 : [200,300],
        2 : [300,450],
        3 : [450, np.inf]
    }
    edges = [0,200,300,450,600]
    #
    sig_p   = ['ttH','ttZ']
    #
    f_dir   = f'{sys.path[1]}/Higgs-Combine-Tool/'
    fit_f   = 'EFT_Parameterization_v3.npy' 
    fit     = np.load(f_dir+fit_f, allow_pickle=True)

    wc_ranges = { 
        # most recent 2sigma constraints
        # if added more points, add after and not inbetween 2 sigma contraints
        'ctp'   : [-10.42,40.24, 9.251, 15, 20, 30, 35],#[-3.97,5.04],
        'cpQM'  : [-9.80,13.22],#[-5.95,11.06],
        'cpQ3'  : [-5.72,5.66],#[-4.913,4.883],
        'cpt'   : [-17.54,12.56],#[-13.50,7.28],
        'cptb'  : [-14.93,14.93],#[-11.47,11.57],
        'ctW'   : [-1.53,1.47],#[-1.19,1.20],
        'ctZ'   : [-1.54,1.57],#[-1.28,1.22],
        'cbW'   : [-6.72,6.72],#[-5.26,5.36],
        'ctG'   : [-0.47,0.36],#[-0.37,0.28],
        'cQei'  : [-200,200],
        'cQl3i' : [-200,200],
        'cQlMi' : [-200,200],
        'ctei'  : [-200,200],
        'ctlSi' : [-200,200],
        'ctlTi' : [-200,200],
        'ctli'  : [-200,200]
      }

    
    @save_pdf('wc_single_effect.pdf') 
    def run_singles(self, years=cfg.Years):
        for y in years:
            for wc, r in self.wc_ranges.items():
                dummy_key = f'y{y}_Zhpt1'
                self.worker(wc, r, dummy_key)
                self.plot_slopes(wc, r, dummy_key)
                
    
    def worker(self, wc, r, d_k):
        # construct pt hist out of fit entries
        pqr = [(wc, wc),('sm',wc),('sm','sm')]
        tth_l, tth_h, *tth_a = [[self.get_eftnorm_v1(k=(f'ttH{i}',d_k),pqr=pqr, v=r[j]) for i in self.pt_bins] for j in range(len(r))]
        ttz_l, ttz_h, *ttz_a = [[self.get_eftnorm_v1(k=(f'ttZ{i}',d_k),pqr=pqr, v=r[j]) for i in self.pt_bins] for j in range(len(r))]
        fig, axs = self.initPlt()

        for lh, _r in zip([ttz_l,ttz_h]+ttz_a, r):
            axs[0].step(x=self.edges,
                        y=np.append(lh,0),
                        where='post', label=f'{wc} {_r}')
        #
        axs[0].set_ylabel(r'(d$\sigma^{EFT}$$/$d${P_T}$) $/$ (d$\sigma^{SM}$$/$d${P_T}$) ')
        self.endPlt(axs[0],'ttZ')
        for lh, _r in zip([tth_l,tth_h]+tth_a, r):
            axs[1].step(x=self.edges,
                        y=np.append(lh,0),
                        where='post', label=f'{wc} {_r}')
        #
        self.endPlt(axs[1],'ttH')
        fig.suptitle(f'{wc} 2$\sigma$ Effects')
        #

    def plot_slopes(self,  wc, r, d_k):
        # now make linear vs quad impacts per wc, using linear regression ot get the slope of the distribution
        fig, axs = self.initPlt()
        tth_i = np.array([[self.get_eftnorm_v1(k=(f'ttH{i}',d_k),pqr=pqr, v=wc_v) for i in self.pt_bins] for wc_v in np.linspace(r[0],r[1], 20+1)])
        ttz_i = np.array([[self.get_eftnorm_v1(k=(f'ttZ{i}',d_k),pqr=pqr, v=wc_v) for i in self.pt_bins] for wc_v in np.linspace(r[0],r[1], 20+1)])
        x = np.array([i for i in self.pt_bins])

        slopes_tth = ((x*tth_i).mean(axis=1) - x.mean()*tth_i.mean(axis=1)) / ((x**2).mean() - (x.mean())**2)
        slopes_ttz = ((x*ttz_i).mean(axis=1) - x.mean()*ttz_i.mean(axis=1)) / ((x**2).mean() - (x.mean())**2)

        axs[0].scatter(x=np.linspace(r[0],r[1], 20+1), y = slopes_ttz, marker='.')
        axs[0].grid(True)
        axs[0].set_ylabel('slope')
        axs[0].set_xlabel(f'{wc}')
        #
        axs[1].scatter(x=np.linspace(r[0],r[1], 20+1), y = slopes_tth, marker='.')
        axs[1].grid(True)
        axs[1].set_xlabel(f'{wc}')
        #
        fig.suptitle(rf'Lin. regression slope of diff. $p_T$ impact vs. {wc}')
        #exit()

    def initPlt(self):
        fig, axs = plt.subplots(1,2, sharey=True)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.0,
            wspace=0.0
        )
        return fig, axs

    def endPlt(self,ax,title=''):
        ax.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True, label='SM')
        ax.set_xlim(0,600)
        ax.set_xlabel(rf"GEN {re.search(r'[ZH]',title).group()} $p_T$ $(GeV)$")
        if title == 'ttZ': ax.legend() # only do once (might just ancor to outside of "box")
        ax.set_title(title)
        #
        tlabels = [str(i) for i in self.edges[:-1]]+['>450']
        if title == 'ttH': tlabels[0] = ''
        ax.set_xticks(self.edges)
        ax.set_xticklabels(tlabels)



    def get_eftnorm_v1(self, k, pqr, v):
        return self.fit[k][pqr[0]]*v*v + self.fit[k][pqr[1]]*v +self.fit[k][pqr[2]]
        

if __name__ == '__main__':
    _  = main()
    _.run_singles(['2018'])


