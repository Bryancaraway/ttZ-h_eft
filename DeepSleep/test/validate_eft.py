import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
from lib.fun_library import save_pdf
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
rc("figure", max_open_warning=600)
import config.ana_cff as cfg

'''
Short script to plot shape validation plots
'''

class main():
    

    sig_p = {'TTZ': ['Zbb','Zqq','Zllnunu'],
             'TTH': ['Hbb','Hnonbb']}
    
    f_dir = f'{sys.path[1]}/files/year/mc_files/'
    years = ['2017','2018'] # for now

    @save_pdf('eft_validation.pdf')
    def run(self):
        for y in self.years:
            for k,v in self.sig_p.items():
                self.worker(k,v,self.f_dir.replace('year',y), y)

    def worker(self, p, sub_p, i_dir, y):
        nom = pd.read_pickle(f'{i_dir}{p}_val.pkl')
        eft = pd.read_pickle(f'{i_dir}{p}_EFT_val.pkl')
        #

        for sp in sub_p:
            fig, ax = plt.subplots()
            cut = (lambda df: (df['NN']>= 0) & (df[sp] == True))
         
            self.plots_with_stats(nom['genZHpt'][cut(nom)].clip(0,500),
                                  nom['genWeight'][cut(nom)]/sum(nom['genWeight'][cut(nom)]),
                                  f'Nom. {sp}',
                                  'orange', ax)
            self.plots_with_stats(eft['genZHpt'][cut(eft)].clip(0,500), 
                                  eft['EFT183'][cut(eft)]/sum(eft['EFT183'][cut(eft)]),
                                  f'EFT {sp}',
                                  'blue', ax)

            ax.legend()
            ax.set_title(f'{p} EFT vs Nom {y}')
            #
            fig, ax = plt.subplots()
            cut = (lambda df: (df['NN']>= 0) & (df[sp] == True))
         
            self.plots_with_stats(nom['genZHpt'][cut(nom)].clip(0,500),
                                  nom['genWeight'][cut(nom)]/sum(nom['genWeight'][cut(nom)]),
                                  f'Nom. {sp}',
                                  'orange', ax)
            self.plots_with_stats(eft['genZHpt'][(eft['EFT183'] < 100) & cut(eft)].clip(0,500), 
                                  eft['EFT183'][(eft['EFT183'] < 100) & cut(eft)]/sum(eft['EFT183'][(eft['EFT183'] < 100) & cut(eft)]),
                                  f'EFT (fix) {sp}',
                                  'blue', ax)
            ax.legend()
            ax.set_title(f'{p} EFT (fix) vs Nom {y}')
            #plt.show() 
            #plt.close('all')
            fig, ax = plt.subplots()
            ax.hist(eft['EFT183'][cut(eft)].clip(0,1000), 
                    #bins=np.quantile(eft['EFT183'][cut(eft)].clip(0,10),np.linspace(0,1,51)), 
                    bins=np.logspace(np.log10(10**(-8)),np.log10(1000), 100),
                    histtype='step',
                    range=(0,1000))
            ax.set_title(f'{p} {sp} EFT183 weight {y}')
            ax.set_yscale('log')
            ax.set_xscale('log')
            #
            #plt.close('all')
            fig, ax = plt.subplots()
            ax.hist((eft['EFT183'][cut(eft)]/sum(eft['EFT183'][cut(eft)])).clip(0,1), 
                    bins=np.logspace(np.log10(10**(-8)),np.log10(1), 100),
                    histtype='step',
                    range=(0,1))
            ax.set_title(f'{p} {sp} EFT183 weight / sum(EFT183) {y}')
            ax.set_yscale('log')
            ax.set_xscale('log')
            #plt.show()
            
            
            
    def plots_with_stats(self, i , w, l, c, ax, b=10, r=(0,500)):
        sumw,  edges = np.histogram(i, weights=w,             bins=b, range=r)
        sumw2, _ = np.histogram(i, weights=np.power(w,2), bins=b, range=r)
        ax.step(    x=edges, 
                     y=np.append(sumw,0), where='post',    color=c, label=l)
        ax.errorbar(x=(edges[1:]+edges[:-1])/2, 
                     y=sumw, yerr=np.sqrt(sumw2), fmt='.', color=c)
        ax.set_xlim(r)
        ax.set_ylim(0,0.4)
        

if __name__ == '__main__':
    _ = main()
    _.run()
