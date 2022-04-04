import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
import numpy as np
np.random.seed(1)
import re
import uproot
from scipy.stats import norm
import config.ana_cff as cfg
from lib.fun_library import t2Run, save_pdf, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

bkgs = ['TTBar','tt_B']
#channels=['Zhpt1','Zhpt2','Zhpt2'] 
channels=['Zhpt2','Zhpt2'] 
dcdir = f'{sys.path[1]}/Higgs-Combine-Tool/'

def do_nui_variations(nui='toppt'):
    nui_hists = {y: get_hists(nui,y) for y in cfg.Years}
    plot_var_by_year(nui_hists, nui)

def get_hists(n,y):
    with uproot.open(dcdir+f'datacard_{y}.root') as roo:
        get_full_hist = (lambda tag : np.sum([np.concatenate([np.array(roo[f'{ch}_{bkg}{tag}'].values)  
                                                              for ch in channels]) for bkg in bkgs], axis=0))
        nom = get_full_hist("")
        up  = get_full_hist(f"_{n}Up")
        dn  = get_full_hist(f"_{n}Down")
        return [up/nom, dn/nom]
    
def plot_var_by_year(hist_dict, nui):
    fig, ax = plt.subplots()
    for c,y in zip(['red','blue','orange'],hist_dict):
        ax.step(x=np.arange(len(hist_dict[y][0])+1), y=np.append(hist_dict[y][0],hist_dict[y][0][-1]),
                where='post', c=f'tab:{c}')
        ax.step(x=np.arange(len(hist_dict[y][0])+1), y=np.append(hist_dict[y][1],hist_dict[y][1][-1]),
                where='post', c=f'tab:{c}', label=f'{nui}_{y}')
    ax.legend()
    ax.set_xlim(0,len(hist_dict['2016'][0]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True,which='both',axis='x')
    #ax.grid(True,axis='y')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    import_mpl_settings()
    do_nui_variations()
