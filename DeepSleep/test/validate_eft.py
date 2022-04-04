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
    

    sig_p = {
        #'TTZ': ['Zbb','Zqq','Zllnunu'],
        #'TTH': ['Hbb','Hnonbb'],
        'TTbb': ['tt_B'],
    }
    sig_k = {'TTZ': 'genZHpt',
             'TTH': 'genZHpt',
             'TTbb': 'ttbb_genbb_pt'}
             #'TTbb': 'Zh_pt'}
    sig_sm = {'TTH': 'ttH',
              'TTZ': 'ttZ',
              'TTbb':'ttbb'}
    
    f_dir = f'{sys.path[1]}/files/year/mc_files/'
    years = ['2016','2017','2018'] # for now

    @save_pdf('eft_validation_ttbb.pdf')
    def run(self):
        for y in ['2016','2017','2018']:
            for k,v in self.sig_p.items():
                self.worker(k,v,self.f_dir.replace('year',y), y)

    def worker(self, p, sub_p, i_dir, y):
        nom    = pd.read_pickle(f'{i_dir}{self.sig_sm[p]}_val.pkl')
        eft    = pd.read_pickle(f'{i_dir}{p}_EFT_val.pkl')
        #eftjet = pd.read_pickle(f'{i_dir}{p}jet_EFT_val.pkl')
        #
        bins = [0,200,300,450,600]
        r = (bins[0],bins[-1])
        #bins = 50 # for ttbb jet study
        #r = (0,250) # for ttbb jet study

        #self.sig_k = {'TTbb':'Zh_pt'}
        #self.sig_k = {'TTbb':'Zh_M'}
        #self.sig_k = {'TTbb':'n_ak4jets'}
        #self.sig_k = {'TTbb':cfg.nn}
        #self.sig_k = {'TTbb':'ttbb_genbb_invm'}
        for sp in sub_p:
            fig, ax = plt.subplots()
            cut = (lambda df: (df[cfg.nn]>= 0) & (df[sp] == True))
            print(f'Raw counts: \nCentral for {p},{y} with {sp}, passing NN cuts: {len(nom[cut(nom)])}\nPrivate for {p},{y} with {sp}, passing NN cuts: {len(eft[cut(eft)])}')
            self.plots_with_stats(nom[self.sig_k[p]][cut(nom)].clip(0,600),
                                  nom['genWeight'][cut(nom)]/sum(nom['genWeight'][cut(nom)]),
                                  f'Cen. {sp}',
                                  'orange', ax,
                                  b=bins,
                                  r=r
                              )
            self.plots_with_stats(eft[self.sig_k[p]][cut(eft)].clip(0,600), 
                                  eft['EFT183'][cut(eft)]/sum(eft['EFT183'][cut(eft)]),
                                  f'Priv. {sp}',
                                  'blue', ax,
                                  b=bins,
                                  r=r
                              )
            #self.plots_with_stats(eftjet[self.sig_k
            #                      eftjet['EFT183'][cut(eftjet)]/sum(eftjet['EFT183'][cut(eftjet)]),
            #                      f'Priv. {sp}+extra parton',
            #                      'red', ax,
            #                      b=bins,
            #                      r=r
            #                  )

            #ax.set_yscale('log')
            ax.legend()

            #ax.set_xlabel(r'GEN bb mass (GeV)')
            ax.set_xlabel(self.sig_k[p])
            
            #ax.set_xlabel(r'Z/H cand. $p_{T}$ (GeV)')
            #ax.set_xlabel(r'Z/H cand. $m_{SD}$ (GeV)')
            #ax.set_xlabel(r'jet multiplicity')
            #ax.set_xlabel(r'MVA score')
            ax.set_ylabel('Fraction of total / bin')
            ax.set_title(f'{p} Private vs Central {y}')
            #plt.show()
            #exit()
            #
            fig, ax = plt.subplots()
            cut = (lambda df: (df[cfg.nn]>= 0) & (df[sp] == True))
         
            self.plots_with_stats(nom[self.sig_k[p]][cut(nom)].clip(0,600),
                                  nom['genWeight'][cut(nom)]/sum(nom['genWeight'][cut(nom)]),
                                  f'Cen. {sp}',
                                  'orange', ax,
                                  b=100
                              )
            self.plots_with_stats(eft[self.sig_k[p]][(eft['EFT183'] < 100) & cut(eft)].clip(0,600), 
                                  eft['EFT183'][(eft['EFT183'] < 100) & cut(eft)]/sum(eft['EFT183'][(eft['EFT183'] < 100) & cut(eft)]),
                                  f'Priv. (fix) {sp}',
                                  'blue', ax,
                                  b=100
                              )
            ax.legend()
            #ax.set_xlabel(r'GEN bb $p_{T}$ (GeV)')
            ax.set_xlabel(self.sig_k[p])
            ax.set_ylabel('Fraction of total / bin')
            ax.set_title(f'{p} Private (weight fix) vs Central {y}')
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
            ax.set_xlabel("EFTrwgt 'SM'")
            ax.set_ylabel('Count / bin')
            #
            #plt.close('all')
            fig, ax = plt.subplots()
            ax.hist((eft['EFT183'][cut(eft)]/sum(eft['EFT183'][cut(eft)])).clip(0,1), 
                    bins=np.logspace(np.log10(10**(-8)),np.log10(1), 100),
                    histtype='step',
                    range=(0,1))
            ax.set_title(f'{p} {sp} EFT183 weight / sum(EFT183) {y}')
            ax.set_xlabel(r"EFTrwgt 'SM' / $\Sigma$ EFTrwgt 'SM'")
            ax.set_ylabel('Count / bin')
            ax.set_yscale('log')
            ax.set_xscale('log')
            #plt.show()
            
            
            
    def plots_with_stats(self, i , w, l, c, ax, b=10, r=(0,600)):
        sumw,  edges = np.histogram(i, weights=w,             bins=b, range=r)
        sumw2, _ = np.histogram(i, weights=np.power(w,2), bins=b, range=r)
        ax.step(    x=edges, 
                     y=np.append(sumw,0), where='post',    color=c, label=l)
        ax.errorbar(x=(edges[1:]+edges[:-1])/2, 
                     y=sumw, yerr=np.sqrt(sumw2), fmt='.', color=c)
        ax.set_xlim(r)
        #ax.set_ylim(0,0.75)
        ax.set_yscale('log')
        if type(b) is list:
            ax.set_xticks(b)
            ax.set_xticklabels([str(i) for i in b[:-1]]+['>450'])


if __name__ == '__main__':
    _ = main()
    _.run()
