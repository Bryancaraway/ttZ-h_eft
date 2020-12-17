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
from matplotlib.ticker import AutoMinorLocator
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg

'''
Script that demonstrates how wc impacts
our expected signal in terms of pt distribution
'''



class main():
    fit_f   = 'EFT_Parameterization_test.npy'
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
    samples=['ttbb','ttjets','ttZ','ttH']
    df = pd.read_pickle(fit_f)

    @save_pdf('ttbb_wc_single_effect.pdf')
    def run_singles(self, years=['2018']):
        norm_change = {}
        for y in years:
            df = self.df[y]['ttH']
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
        plt.show()
        exit()
    
    #@save_pdf('eft_bkgsig_comparison_nn09.pdf')
    def run_singles_test(self, samples, cut, title, year='2018' ):
        for wc, r in self.wc_ranges.items():
            fig,ax = self.initPlot()
            dott2b = True
            for s in samples:
                kinem = 'Zh_pt'
                
                #sam_text = f'{s}'
                df = self.df[year][s]
                slabel = s
                df = cut(df)
                if s == 'ttbb':
                    if dott2b:
                        df = df[df['process'] == 'tt_2b']
                        slabel = 'tt_2b'
                        dott2b=False
                    else: 
                        df = df[df['process'] == 'tt_bb']
                elif s == 'ttjets':
                    df = df[df['process'] == 'TTBar']
            
                nsm, b= np.histogram(
                    df[kinem].clip(0,590),
                    bins=[200,300,450,600],
                    weights=df['SM']
                )
                #
                tot_sm, _  = np.histogram(
                    df[kinem],
                    bins=[200,300,450,600],
                    weights=df['SM']*df['weight']
                )
                print(slabel,len(df))
                #
                tot_eft = {0:[],1:[]}
                #
                for i,v in enumerate(r):
                    eftw = self.get_eftw(df, wc, v)
                    neft, _ = np.histogram(
                        df[kinem].clip(0,590),
                        bins=[200,300,450,600],
                        weights=eftw
                    )
                    tot_eft[i] = (neft/nsm)*tot_sm
                    ax.step(
                        x = b,
                        y = np.append(neft/nsm,0),
                        where='post', label=f'{slabel}:{v}'
                    )
                #for i,xy in enumerate(zip((b[1:]+b[:-1])/2,neft/nsm)):
                #ax.text(x=xy[0],y=xy[1], s=f'{slabel} {tot_sm[i]:.1f}:{tot_eft[0][i]:.1f}:{tot_eft[1][i]:.1f}', fontsize=8)
            ax.legend(ncol=2, fontsize='x-small')
            ax.set_xlim(200,600)
            ax.set_xlabel(' reco. Z/H_pt (GeV)')
            ax.set_ylabel(r'(d$\sigma^{EFT}$$/$d${P_T}$) $/$ (d$\sigma^{SM}$$/$d${P_T}$) ')
            #fig.suptitle(f'EFT impact, NN > 0.9, {wc}:{r}')
            fig.suptitle(title.format(wc=wc, r=r))
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            #plt.yscale('log')
            #plt.show()
            #exit()
            
    def initPlot(self):
        fig, ax = plt.subplots(1,1, sharey=True)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.0,
            wspace=0.0
        )
        return fig,ax

    def get_eftw(self, df, wc, v):
        p = df[f'{wc}_{wc}']*v*v
        q = df[f'{wc}']*v
        r = df['SM']
        _eft_w = p + q + r
        return _eft_w

    def calc_norm(self,df,wc,v):
        p = sum(df[f'{wc}_{wc}']*v*v)
        q = sum(df[f'{wc}']*v)
        r = sum(df['SM'])
        return p/r + q/r + r/r
    
    def getMCStat_info(self,n, b, h, w):
        #x = (self.bins[1:]+self.bins[:-1])/2
        x = (b[1:]+b[:-1])/2
        xerr = (b[1:]-b[:-1])/2
        yerr = np.histogram(h, bins=b, weights=np.power(w,2))[0]
        yerr = np.sqrt(yerr)
        y = n
        return x,xerr,y,yerr
    
if __name__=='__main__':
    _ = main()
    #_.run_singles()
    cut_nntight = (lambda x : x[(x['NN'] > .9)])
    cut_nnloose = (lambda x : x[(x['NN'] > 0)])
    #_.run_singles_test(['ttH','ttZ','ttbb','ttbb'],cut)
    #@save_pdf('eft_bkgsig_comparison_nn09.pdf')
    save_pdf('eft_bkgsig_comparison_nn09.pdf')(_.run_singles_test)(['ttH','ttZ','ttbb','ttbb'],cut_nntight, 'EFT impact, NN > 0.9, {wc}:{r}')
    save_pdf('eft_bkgsig_comparison_nn00.pdf')(_.run_singles_test)(['ttH','ttZ','ttbb','ttbb','ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')
    save_pdf('eft_bkg_comparison_nn09.pdf')(_.run_singles_test)(['ttbb','ttbb'],cut_nntight, 'EFT impact, NN > 0.9, {wc}:{r}')
    save_pdf('eft_bkg_comparison_nn00.pdf')(_.run_singles_test)(['ttbb','ttbb','ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')
    
