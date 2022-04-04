import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import re
import seaborn as sns
import config.ana_cff as cfg
#from lib.fun_library import save_pdf, getLaLabel
from lib.fun_library import save_pdf, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel
from qcDatacard import tbins_map
#from post_fit import PostFit
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
#from matplotlib.collections import PatchCollection
#from matplotlib.patches import Patch, Rectangle
#from matplotlib import transforms
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)

p_dict = {'ttbb':'tt_B',
          'TTBar':'TTBar'}

@save_pdf('ttpt_investigation.pdf')
def main():
    # get top data
    ttbb = get_val_info('ttbb')
    ttbar  = get_val_info('TTBar')
    data = pd.concat([ttbb,ttbar], axis='rows')
    # define new mass bins
    data['mass_bin'] = ''
    data.loc[(data['Zh_M']>75) & (data['Zh_M']<105),'mass_bin'] = 'Z mass bin'
    data.loc[(data['Zh_M']>105) & (data['Zh_M']<145),'mass_bin'] = 'H mass bin'
    data = data.loc[(data['mass_bin'].str.contains('mass')) & (data['Zh_pt'] < 600)]
    data['tt_sum_pt'] = data['tt_pt1']+data['tt_pt2']
    data['NN_score'] = 'All events'
    data['NN_score'].loc[data[cfg.nn]>0.8] = 'NN > 0.80'
    # scatter versus Zh_pt, Zh_mass
    #fig, ax = plt.subplots()
    #CMSlabel(fig,ax,altloc=False,opt='Simulation Preliminary', lumi='nl')
    #ax.scatter(x=ttbb['topptWeight'],y=ttbb['Zh_M'],c='tab:orange',label='ttbb')
    #ax.scatter(x=ttbar['topptWeight'],y=ttbar['Zh_M'],c='tab:blue',label='TTbar')
    plot = sns.displot(data,x='Zh_pt',y='topptWeight',hue='mass_bin', kind='kde', levels=5)
    plot.set(title="All")
    plt.tight_layout()
    plot = sns.displot(data[data['NN_score']=='NN > 0.80'],x='Zh_pt',y='topptWeight',hue='mass_bin', kind='kde', levels=5)
    plot.set(title="NN > 0.8")
    plt.tight_layout()
    #fig, ax = plt.subplots()
    plot = sns.displot(data,x='Zh_pt',y='tt_sum_pt',hue='mass_bin', kind='kde', levels=5)
    plot.set(title="All")
    plt.tight_layout()
    plot = sns.displot(data[data['NN_score']=='NN > 0.80'],x='Zh_pt',y='tt_sum_pt',hue='mass_bin', kind='kde', levels=5)
    plot.set(title="NN > 0.8")
    plt.tight_layout()
    #ax.set_xlabel(r'topptWeight')
    #ax.set_ylabel(r'Z/H $m_{\text{SD}}$')
    #ax.legend()

    #plt.show()

def get_val_info(sample, y='2018'):
    s_df = pd.read_pickle(f'{sys.path[1]}/files/{y}/mc_files/{sample}_val.pkl').filter(items=['tt_pt1','tt_pt2','topptWeight','Zh_pt','Zh_M','process', cfg.nn])
    return s_df[(s_df[cfg.nn] > 0.0) & (s_df['process'] == p_dict[sample])]

if __name__ == '__main__':
    #import_mpl_settings(2,disable_sansmath=True)
    main()
