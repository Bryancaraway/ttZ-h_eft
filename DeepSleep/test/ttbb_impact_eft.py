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
    fit_f   = 'EFT_Parameterization_ttbb.npy'
    wc_ranges = {
        'ctW'  :[-2.44, 2.40]  ,
        'ctZ'  :[-2.33, 2.31]  ,
        'ctp'  :[-11.64, 41.58],
        'cpQM' :[-12.33, 16.04],
        'ctG'  :[-0.50, 0.48]  ,
        'cbW'  :[-7.40, 7.40]  ,
        'cpQ3' :[-6.94, 7.20]  ,
        'cptb' :[-16.67, 16.67],
        'cpt'  :[-20.35, 16.53]
        }
    sample='ttbb'
    df = pd.read_pickle(fit_f)

    @save_pdf('ttbb_wc_single_effect.pdf')
    def run_singles(self, years=['2018']):
        norm_change = {}
        for y in years:
            df = self.df[y]['ttbb']
            for wc, r in self.wc_ranges.items():
                norm_change[wc] = [self.calc_norm(df,wc,v) for v in r]
    
        wcs   = np.array([w for w in norm_change])
        norms = np.array([norm_change[w] for w in norm_change])
        ind   = np.array([ [i,i] for i in range(len(norm_change))])
        for i, w in enumerate(norm_change):
            plt.scatter([i,i],norm_change[w],c='k')
            plt.plot([i,i],norm_change[w],c='k')
        #plt.scatter(ind.flatten(), norms.flatten())
        plt.axhline(1,c='r',ls='--')
        plt.fill_between([-10,10],1.21,0.82,color='r',alpha=0.5,label=f'Exp $1\sigma$ run2 limit')
        plt.fill_between([-10,10],1.42,0.64,color='k',alpha=0.25,label=f'Exp $2\sigma$ run2 limit')
        plt.xticks([i for i in range(len(norm_change))], wcs)
        plt.ylabel(r'$\sigma$(ttbb)$^{EFT}$/$\sigma$(ttbb)$^{SM}$')
        plt.xlim(-1,10)
        plt.title('Rate change w/resp. to WC for tt+bb')
        plt.legend()
        #plt.show()
        

    def calc_norm(self,df,wc,v):
        p = sum(df[f'{wc}_{wc}']*v*v)
        q = sum(df[f'{wc}']*v)
        r = sum(df['SM'])
        return p/r + q/r + r/r
    
if __name__=='__main__':
    _ = main()
    _.run_singles()
