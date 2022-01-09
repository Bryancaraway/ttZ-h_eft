import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import matplotlib.pyplot as plt
from matplotlib import rc
rc("figure", max_open_warning=600)
#rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

import json
import numpy as np
import pandas as pd
import config.ana_cff as cfg
from lib.fun_library import t2Run, save_pdf, getZhbbBaseCuts, getZhbbWeight, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel

nn = cfg.nn
year = '2018'

def check_unc(unc_name, sample, cut_func):
    p_df = cut_func( get_sample(sample) )
    p_weight = getZhbbWeight(p_df, year)
    up,down = sum(p_weight*p_df[unc_name+'_Up'])/sum(p_weight), sum(p_weight*p_df[unc_name+'_Down'])/sum(p_weight)
    print(f"{sample} rate variation for {unc_name} | Up: {up:.2f} ({abs(1-up)*100:.2f}%) | Down: {down:.2f} ({abs(1-down)*100:.2f}%)\n")

def get_sample(sample):
    p_file = f'{cfg.master_file_path}/{year}/mc_files/{sample}_val.pkl'
    p_df = pd.read_pickle(p_file)
    return p_df

if __name__ == '__main__':
    #check_unc('ISR', 'TTBar', (lambda x: x[( (x[cfg.nn] >= 0.0) & (x['process'] == 'TTBar') & (x['tt_C'] == True) )])) 
    #check_unc('FSR', 'TTBar', (lambda x: x[( (x[cfg.nn] >= 0.0) & (x['process'] == 'TTBar') & (x['tt_C'] == True) )])) 
    ##
    #check_unc('ISR', 'ttbb', (lambda x: x[( (x[cfg.nn] >= 0.0) & (x['process'] == 'tt_B')  )])) 
    #check_unc('FSR', 'ttbb', (lambda x: x[( (x[cfg.nn] >= 0.0) & (x['process'] == 'tt_B')  )])) 
    # mu_r_Up/Down, mu_f_Up/Down
    check_unc('mu_r', 'ttbb', (lambda x: x[( (x[cfg.nn] >= 0.0) & (x['process'] == 'tt_B') & (x['tt_2b'] == True)  )])) 
    check_unc('mu_f', 'ttbb', (lambda x: x[( (x[cfg.nn] >= 0.0) & (x['process'] == 'tt_B') & (x['tt_2b'] == True)  )])) 
    
