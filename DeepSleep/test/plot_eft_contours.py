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
from matplotlib import rc, lines
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
wc_bins = {'ctp':np.arange(-15, 50, 1),
           'cpt':np.arange(-20, 20, .5),
           'ctZ':np.arange(-3,3, .2),} # cant do -8, 8
norm_bins = {'cpt':np.arange(1,5,.1),
             'ctZ':np.arange(1,5,.1)}
#>>> a.max(), a.min()
#ctW       5.949762   -5.890477
#ctp      61.349454  -16.477814
#cpQM     41.893451  -11.870914
#ctZ       5.918831   -5.982110
#cbW       7.492276   -7.454655
#cpQ3      5.901889  -14.802903
#cpt      19.332344  -35.709867
#cptb     26.947300  -26.654921



p_to_wc = {
    'ttZ':'ctZ', 
    'ttH':'cpt'}
wc_latex = {
    #'cbW'  : r'${c}_{\mathrm{bW}}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'cptb' : r'${c}_{\varphi \mathrm{tb}}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'cpt'  : r'${c}_{\varphi \mathrm{t}}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'ctp'  : r'${c}_{\mathrm{t} \varphi}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'ctZ'  : r'${c}_{\mathrm{tZ}}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'ctW'  : r'${c}_{\mathrm{tW}}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'cpQ3' : r'${c}_{\varphi \mathrm{Q}}^{3}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    #'cpQM' : r'${c}_{\varphi \mathrm{Q}}^{-}\,/\,{\Lambda}^{2} \; [{\mathrm{TeV}}^{-2}]$',
    'ctp'  : r'$\mathsf{c_{t \varphi}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cpQM' : r'$\mathsf{c^{-}_{\varphi Q}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cpQ3' : r'$\mathsf{c^3_{\varphi Q}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cpt'  : r'$\mathsf{c_{\varphi t}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cptb' : r'$\mathsf{c_{\varphi t b}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'ctW'  : r'$\mathsf{c_{tW}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'cbW'  : r'$\mathsf{c_{bW}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
    'ctZ'  : r'$\mathsf{c_{tZ}} \,/\, \Lambda^\mathsf{2} $ \raisebox{0.25ex}{[}$\text{TeV}^\mathsf{-2}$\raisebox{0.25ex}{]}',
}

#@save_pdf("eft_process_contours_nolumi_final.pdf")
@save_pdf("eft_process_contours_nolumi.pdf")
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
        eft_pt_zbin_content, norm_pt_zbin_content = bin_2d_eft_effects(p_df) 
        plot_eft_y_contours(x= (pt_bins[1:]+pt_bins[:-1])/2, 
                           y= wc_bins[p_to_wc[p]], 
                           z= eft_pt_zbin_content,
                           process= p)
        plot_eft_z_contours(x= (pt_bins[1:]+pt_bins[:-1])/2, 
                            y= norm_bins[p_to_wc[p]], 
                            z= norm_pt_zbin_content,
                            process= p)
        
        #
    #
    
def plot_eft_y_contours(x, y, z, process):
    fig, ax = beginPlt()
    #levels = [2,5,10]
    levels = [1.2,1.5,2.0]
    #tslabels = [r'$\frac{\sigma_{\mathrm{EFT}}}{\sigma_{\mathrm{SM}}}=$'+str(i) 
    #            for i in levels]
    tslabels = [r'$\frac{\mathsf{\sigma}_{\text{EFT}}}{\mathsf{\sigma}_{\text{SM}}}\mathsf{=}$'+str(i) 
                for i in levels]
    #triang = tri.Triangulation(x, y) # broken
    #ts = ax.tricontour(x, y , z, levels=levels, colors=['gold','blue','green','magenta']) # broken
    X,Y = np.meshgrid(x,y)
    ts = ax.contour(X, Y , z, levels=levels, colors=['orange','green','blue'], linewidths=0.5)
    ax.axhline(0, color='red', linewidth=.5, snap=True)
    # === works filled contour
    #cntr = ax.contourf(x, y , z, levels=np.arange(0,5+.5,.5))
    #cbar = fig.colorbar(cntr, pad=.05)
    #cbar.ax.set_ylabel(r'$\sigma_{\mathrm{EFT}}\,/\,\sigma_{\mathrm{SM}}$')
    # === end of working filled contour

    #ax.clabel(ts, fmt={l:ls for l,ls in zip(levels, tslabels)}, inline=1, fontsize=8, manual=False)
    # setup labels, handles for legend
    handles = [
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='red'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='orange'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='green'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='blue'),
    ]
    labels = ['1.0' , '1.2' , '1.5', '2.0']
    leg = ax.legend(handles, labels, handlelength=1.0, fontsize=6.0, ncol=len(labels), framealpha=0, 
                    loc='upper right',
                    bbox_to_anchor=(1.01, 1.05) if p_to_wc[process] == 'ctZ' else (.55,.8),
                    #title=r'$\frac{\sigma_{\mathrm{EFT}}}{\sigma_{\mathrm{SM}}}=$' )
                    #title=r'$\sigma_{\mathrm{EFT}}/\sigma_{\mathrm{SM}}=$' )
                    title=r'$\mathsf{\sigma}_{\text{EFT}}/\mathsf{\sigma}_{\text{SM}}\mathsf{=}$' )
    leg._legend_box.align = 'left'
    #
    ax.set_xlim(0,600)
    #for i in range(len(tslabels)):
    #    ts.collections[i].set_label(tslabels[i])
    z_or_h = re.search(r'(Z|H)',process).group()
    #ax.set_xlabel(rf"Simulated ${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{z_or_h}}}}}$ [GeV]")
    #ax.set_xlabel(rf"${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{z_or_h}}}}}$ [GeV]")
    ax.set_xlabel(rf"$\mathsf{{p}}_{{\text{{T}}}}^{{\text{{{z_or_h}}}}} \; \left[\smash{{\text{{GeV}}}}\right]$", usetex=True)
    ax.set_ylabel(wc_latex[p_to_wc[process]])
    plt.xlim(50,600)
    plt.tight_layout()
    #plt.show()
    #exit()

def plot_eft_z_contours(x, y, z, process):
    fig, ax = beginPlt()
    #levels = [2,5,10]
    level_dict = {'ctZ':[0.6,1.0,1.2],
                  'cpt':[5.0,8.0,10.0]}
    levels = level_dict[p_to_wc[process]]
    tlabel = wc_latex[p_to_wc[process]]
    tslabels = [tlabel+r'$=$'+str(i) for i in levels]
    X,Y = np.meshgrid(x,y)
    import scipy.ndimage
    #z = scipy.ndimage.zoom(z,3)
    z = scipy.ndimage.filters.gaussian_filter(z, sigma=.6)
    #X = scipy.ndimage.zoom(X,3)
    #Y = scipy.ndimage.zoom(Y,3)
    ts = ax.contour(X, Y , z, levels=levels, colors=['orange','green','blue'], linewidths=0.5, antialiased=True,)
    #ts = ax.contour(z, levels=levels, colors=['orange','green','blue'], linewidths=0.5, antialiased=True,)
    ax.axhline(1, color='red', linewidth=.5, snap=True)
    # setup labels, handles for legend
    handles = [
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='red'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='orange'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='green'),
        lines.Line2D([],[], linestyle='-', linewidth=.5, color='blue'),
    ]
    #labels = ['1.0' , '1.2' , '1.5', '2.0']
    labels = ["0.0"]+[f'{i:.1f}' for i in levels]
    leg = ax.legend(handles, labels, handlelength=1.0, fontsize=10.0, ncol=len(labels), framealpha=0, 
                    loc='upper left',
                    bbox_to_anchor=(-.02, 1.07),
                    title=tlabel+r'$\mathsf{=}$')

    leg._legend_box.align = 'left'
    #
    ax.set_xlim(0,600)
    ax.set_ylim(0.8,2)
    z_or_h = re.search(r'(Z|H)',process).group()
    #ax.set_xlabel(rf"${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{z_or_h}}}}}$ [GeV]")
    #ax.set_ylabel(r'$\sigma_{\mathrm{EFT}}/\sigma_{\mathrm{SM}}$')
    #ax.set_xlabel(rf"$\mathsf{{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{z_or_h}}}}} \; \left[\mathrm{{GeV}}\right]$")
    #ax.set_xlabel(rf"$\mathsf{{p}}_{{\text{{T}}}}^{{\text{{{z_or_h}}}}} \; \left[\smash{{\text{{GeV}}}}\right]$", usetex=True)
    ax.set_xlabel(rf"$\mathsf{{p}}_{{\text{{T}}}}^{{\text{{{z_or_h}}}}} $ \raisebox{{0.25ex}}{{[}}$\text{{GeV}}$\raisebox{{0.25ex}}{{]}}")
    ax.set_ylabel(r'$\mathsf{\sigma}_{\text{EFT}}/\mathsf{\sigma}_{\text{SM}}$', usetex=True)
    plt.xlim(50,600)
    plt.tight_layout()
    #plt.show()
    #exit()

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
    #pt_vs_eft = np.array([[calc_norm(wc_v, df=df_group, wc=p_wc) for wc_v in wc_bins[p_wc] ] for df_group in df_groups])
    pt_vs_eft  = np.array([[calc_norm(wc_v, df=df_group, wc=p_wc) for wc_v in wc_bins[p_wc] ] for df_group in df_groups])
    pt_vs_norm = np.array([[calc_wc_from_norm(norm_v, df=df_group, wc=p_wc) for norm_v in norm_bins[p_wc] ] for df_group in df_groups])
    #print(pt_vs_eft.shape)
    return pt_vs_eft.T, pt_vs_norm.T


def beginPlt():
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.88,wspace=0.0,hspace=0.0)
    CMSlabel(fig,ax,altloc=False,opt='Simulation Preliminary', lumi='nl')
    #CMSlabel(fig,ax,altloc=False,opt='Simulation', lumi='nl')
    return fig, ax

class getBeta(TestEFTFitParams):
    def __init__(self, sample):
        self.aux_df = {sample : pd.read_pickle(f'{self.aux_dir}/{self.aux_dict[sample]}')} 

def calc_impact(v,df=None,wc=None):
    p = sum(df[f'{wc}_{wc}']*v*v)
    q = sum(df[f'{wc}']*v)
    r = sum(df['SM'])
    return r/p + r/q + r/r

def calc_norm(v,df=None,wc=None):
    p = sum(df[f'{wc}_{wc}']*v*v)
    q = sum(df[f'{wc}']*v)
    r = sum(df['SM'])
    return p/r + q/r + r/r

def calc_wc_from_norm(norm, df=None, wc=None):
    a = sum(df[f'{wc}_{wc}'])/sum(df['SM'])
    b = sum(df[f'{wc}'])/sum(df['SM'])
    c = sum(df['SM'])/sum(df['SM']) - norm
    return (-1*b + np.sqrt( np.power(b,2) - 4*a*c ))/(2*a)

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
    import_mpl_settings(disable_sansmath=True)
    main()
