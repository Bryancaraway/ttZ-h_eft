import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import uproot
import re
from functools import partial
from multiprocessing import Pool
import concurrent.futures
executor = concurrent.futures.ThreadPoolExecutor()
from lib.fun_library import save_pdf, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
from modules.eftParam import TestEFTFitParams

EFT_param = 'EFT_Parameterization_test.npy'
dummy_y = '2018'

i_file_dict = {
    'ttZ':"files/ttz_eft_2018.txt",
    'ttH':"files/tth_eft_2018.txt"
}

pt_name = 'genZHpt'
pt_bins = np.linspace(0,1000,15)
wc_bins = {'ctp':np.arange(-30, 60, 1),
           'ctZ':np.arange(-8,8, .1),}
p_to_wc = {'ttZ':'ctZ', 'ttH':'ctp'}
wc_latex = {
    'cbW'  : r'${c}_{\mathrm{bW}}\,/\,{\Lambda}^{2}$',
    'cptb' : r'${c}_{\phi \mathrm{tb}}\,/\,{\Lambda}^{2}$',
    'cpt'  : r'${c}_{\phi \mathrm{t}}\,/\,{\Lambda}^{2}$',
    'ctp'  : r'${c}_{\mathrm{t} \phi}\,/\,{\Lambda}^{2}$',
    'ctZ'  : r'${c}_{\mathrm{tZ}}\,/\,{\Lambda}^{2}$',
    'ctW'  : r'${c}_{\mathrm{tW}}\,/\,{\Lambda}^{2}$',
    'cpQ3' : r'${c}_{\phi \mathrm{Q}}^{3}\,/\,{\Lambda}^{2}$',
    'cpQM' : r'${c}_{\phi \mathrm{Q}}^{-}\,/\,{\Lambda}^{2}$',
}

@save_pdf("eft_process_contours.pdf")
def main():
    #eft_dict = np.load(EFT_param, allow_pickle=True)
    #
    for p in p_to_wc:
        eft_df = pd.DataFrame()
        with open(i_file_dict[p]) as ifile:
            eft_df = __worker(ifile.readlines())
        eft_df['process'] = p
        #p_df = pd.concat([eft_dict[y][p] for y in ['2016','2017','2018']], axis='rows', ignore_index=True) 
        p_df = getBeta(p).calcBeta(eft_df,p)  
        del eft_df
        eft_pt_zbin_content = bin_2d_eft_effects(p_df) 
        plot_eft_contours(x= (pt_bins[1:]+pt_bins[:-1])/2, 
                           y= wc_bins[p_to_wc[p]], 
                           z= eft_pt_zbin_content,
                           process= p)
        #
    #
    
def plot_eft_contours(x, y, z, process):
    fig, ax = beginPlt()
    levels = [2,3,5]
    tslabels = [r'$\frac{\sigma_{\mathrm{EFT}}}{\sigma_{\mathrm{SM}}}=$'+str(i) for i in levels]
    #triang = tri.Triangulation(x, y)
    #ts = ax.tricontour(x, y , z, levels=levels, colors=['gold','blue','green','magenta'])
    #X,Y = np.meshgrid(x,y)
    #ts = ax.contour(X, Y , z, levels=levels, colors=['gold','blue','green'], linewidths=0.5)
    cntr = ax.contourf(x, y , z, levels=np.arange(0,5+.5,.5))
    cbar = fig.colorbar(cntr, pad=.05)
    cbar.ax.set_ylabel(r'$\sigma_{\mathrm{EFT}}\,/\,\sigma_{\mathrm{SM}}$')
    
    #ax.clabel(ts, fmt={l:ls for l,ls in zip(levels, tslabels)}, inline=1, fontsize=8, manual=False)
    #ax.set_xlim(0,600)
    #for i in range(len(tslabels)):
    #    ts.collections[i].set_label(tslabels[i])
    z_or_h = re.search(r'(Z|H)',process).group()
    ax.set_xlabel(rf"Simulated ${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{z_or_h}}}}}$ [GeV]")
    ax.set_ylabel(wc_latex[p_to_wc[process]])
    plt.xlim(50,600)
    plt.tight_layout()

def bin_2d_eft_effects(df):
    # create new column with pt bin labels
    df['pt_eft_bin'] = pd.cut(
        df[pt_name].clip(pt_bins[0],pt_bins[-1]), 
        bins=pt_bins,
        labels=np.arange(len(pt_bins)-1)
    )
    # get pt x eft norm matrix
    df_groups = [df[df['pt_eft_bin']==i_pt] for i_pt in sorted(df['pt_eft_bin'].unique())]
    p_wc = p_to_wc[df['process'].iloc[0]]
    pt_vs_eft = np.array([[calc_norm(wc_v, df=df_group, wc=p_wc) for wc_v in wc_bins[p_wc] ] for df_group in df_groups])
    return pt_vs_eft.T

def beginPlt():
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.88,wspace=0.0,hspace=0.0)
    CMSlabel(fig,ax)
    return fig, ax

class getBeta(TestEFTFitParams):
    def __init__(self, sample):
        self.aux_df = {sample : pd.read_pickle(f'{self.aux_dir}/{self.aux_dict[sample]}')} 

def calc_norm(v,df=None,wc=None):
    p = sum(df[f'{wc}_{wc}']*v*v)
    q = sum(df[f'{wc}']*v)
    r = sum(df['SM'])
    return p/r + q/r + r/r

def __worker(lines):
    pool = Pool(8)
    results = pool.map(get_eft_df,[line.strip('\n').replace('root://kodiak-se.baylor.edu//','/cms/data/') for line in lines])
    pool.close()
    _out_df = pd.concat(results, axis='rows', ignore_index=True)
    return _out_df

def get_eft_df(roofile):
    with uproot.open(roofile) as roo:
        t = roo['Events']
        # get eft weights
        eft_reweight = t.array('LHEReweightingWeight', executor=executor)
        # get gen Z/H pt
        g_ids = t.array('GenPart_pdgId', executor=executor)
        g_st = t.array('GenPart_status', executor=executor)
        g_pt = t.array('GenPart_pt', executor=executor)
        #
        isZH = ((g_ids == 25) | (g_ids == 23))
        #
        _df = pd.DataFrame()
        _df['genZHpt'] = g_pt[isZH][:,0].flatten()
        for i in range(184):
            _df[f'EFT{i}'] = eft_reweight[:,i]
        return _df

if __name__ == '__main__':
    import_mpl_settings()
    main()
