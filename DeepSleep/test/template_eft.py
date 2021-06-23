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
from lib.datacard_shapes import DataCardShapes
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg


'''
WC impact per template bin
'''
pt_bins = [0,200,300,450]
sdM_bins  = cfg.sdm_bins
nn = cfg.nn

hist2d = DataCardShapes(pt_bins,sdM_bins,n_NN_bins=10, isblind=True) 

def get_sumw_sumw2(df, weights, year):
    ret_x_y = (lambda df: (df[nn].to_numpy(), df['Zh_M'].to_numpy()))
    sumw = np.array([
        hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
                             weights=weights[df['pt_bin']==i_bin])[0] 
        for i_bin in range(1,len(pt_bins))
    ])
    sumw2 = np.array([
        hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
                             weights=np.power(weights[df['pt_bin']==i_bin],2))[0]
        for i_bin in range(1,len(pt_bins))
    ])
    return sumw, sumw2
#

class Template_EFT:
    fit_f   = 'EFT_Parameterization_test.npy'    

    def __init__(self):
        self.df = pd.read_pickle(self.fit_f)
        df_cut = (lambda _df : _df[_df[nn]>0.0])
        # add 'pt_bin' to df s 
        for y in self.df:
            for s in self.df[y]:
                if len(self.df[y][s]) > 0:
                    self.df[y][s]['pt_bin'] = pd.cut(self.df[y][s]['Zh_pt'].clip(pt_bins[1],pt_bins[-1]+1), 
                                                     bins=pt_bins[1:]+[500],
                                                     labels=[i_bin for i_bin in range(1,len(pt_bins))])
                    self.df[y][s] = df_cut(self.df[y][s])

    def plot_templates(self,procs,years):
        fig,ax = self.initPlot()
        for proc in procs:
            for year in years:
                sumw, sumw2 = get_sumw_sumw2(self.df[year][proc], np.ones(len(self.df[year][proc])), year)
                print(sumw)
                # merge last 2 mass bins of pt bin 1
                sumw[0,:,-2] = sumw[0,:,-2] + sumw[0,:,-1]
                #
                sumw = np.append(sumw[0,:,:-1].flatten(), sumw[1:,:,:].flatten())
                print(sumw)
                print(len(sumw))
                bins = np.arange(len(sumw)+1)
                ax.step(
                    x = bins,
                    y = np.append(sumw,0),
                    where= 'post', label=f'{proc}_{year}')
                ax.errorbar(**self.errbar_kwargs(sumw,np.sqrt(sumw),bins))
                ax.set_yscale('log')
                plt.show()
                exit()


    def errbar_kwargs(self, y, yerr, bins):
        kwargs = {'x':(bins[1:]+bins[:-1])/2, 'y':y, 'yerr':yerr, 'fmt':'.'}
        return kwargs

    #def get_sumw_sumw2(self,df_,w_,bins):
    #    sumw, _ = np.histogram(df_, bins=bins, weights=w_)
    #    sumw2, _ = np.histogram(df_, bins=bins, weights=np.power(w_,2))
    #    return sumw, sumw2
            
    def initPlot(self):
        fig, ax = plt.subplots(1,1)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.0,
            wspace=0.0
        )
        return fig,ax

    def endPlot(self,ax,title,kinem, bins):
        ax.legend(ncol=2, fontsize='x-small')
        ax.set_xlim(bins[0],bins[-1])
        ax.set_ylim(0.4)
        klabel[kinem](ax)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

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


def main():
    template_eft = Template_EFT()
    template_eft.plot_templates(['ttH'],['2018'])
    pass

if __name__ == '__main__':
    main()
