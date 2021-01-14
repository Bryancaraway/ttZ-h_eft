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
import concurrent.futures
executor = concurrent.futures.ThreadPoolExecutor()
from lib.fun_library import save_pdf
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg
from modules.eftParam import TestEFTFitParams



def main(input_file):
    eft_df = pd.DataFrame()
    with open(input_file) as ifile:
        for line in ifile.readlines():
            print(line)
            eft_df = pd.concat([eft_df,get_eft_df(line.strip('\n'))], axis = 'rows', ignore_index=True)
            #
        #
    # get parameterization
    df = getBeta('ttjets').calcBeta(eft_df,'ttjets')
    plot_overal_norm(df, 'ttjets')
    
    
def plot_overal_norm(df, sample):
    wc_ranges = {'ctG'  :np.arange(-1,1,.01) }
    norm_change = {}
    for wc, r in wc_ranges.items():
        norm_change[wc] = [calc_norm(df,wc,v) for v in r]
    
    for w in norm_change:
        fig, ax = plt.subplots()
        ax.plot(wc_ranges[w],norm_change[w])
        ax.axhline(1,c='r',ls='--')
        #plt.xticks([i for i in range(len(norm_change))], norm_change.keys())
        ax.set_xlabel(f'WC {w} ')
        ax.set_ylabel(r'$\sigma$(ttjets)$^{EFT}$/$\sigma$(ttjets)$^{SM}$')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        fig.suptitle(f'Inc. Rate change w/resp. to {w} for tt+jets')
        #plt.xlim(-1,1)
        #plt.legend()
        plt.show()

def get_eft_df(roofile):
    with uproot.open(roofile) as roo:
        t = roo['Events']
        eft_reweight = t.array('LHEReweightingWeight', executor=executor)
        _df = pd.DataFrame()
        for i in range(184):
            _df[f'EFT{i}'] = eft_reweight[:,i]
        return _df
                

class getBeta(TestEFTFitParams):
    def __init__(self, sample):
        self.aux_df = {sample : pd.read_pickle(f'{self.aux_dir}/{self.aux_dict[sample]}')} 

def calc_norm(df,wc,v):
    p = sum(df[f'{wc}_{wc}']*v*v)
    q = sum(df[f'{wc}']*v)
    r = sum(df['SM'])
    return p/r + q/r + r/r


if __name__ == '__main__':
    main("ttjets_eft.txt")
