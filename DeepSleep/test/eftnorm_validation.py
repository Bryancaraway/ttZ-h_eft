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
from lib.fun_library import save_pdf
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg
from modules.eftParam import TestEFTFitParams

single_points = {
    'ttZ': {
        'SM': [[0,1]], 'cbW':[[-7,.6903/.6320],[7,.7009/.6320]], 'cpQ3': [[ -7,.6555/.6320], [ 7,.6909/.6320] ], 'cpQM': [[ -12,1.778/.6320], [ 16,.2684/.6320] ],
        'cpt': [[ -20,.8176/.6320], [ 15,1.852/.6320] ], 'cptb': [[ -15,.6567/.6320], [ 15,.6553/.6320] ], 'ctW': [[ -2, .6796/.6320], [ 2, .7014/.6320] ], 'ctZ': [[ -2,.8729/.6320], [ 2,.8822/.6320] ],
        'ctp': [[ -11,.6304/.6320], [ 40,.6397/.6320] ]
    },
    'ttH' : {
        'SM': [[0,1]], 'cbW':[[-7,.4911/.4461],[7,.4892/.4461]], 'cpQ3': [[ -7,.5008/.4461], [ 7,.5103/.4461] ], 'cpQM': [[ -12, .4954/.4461], [ 16,.5178/.4461] ],
        'cpt': [[ -20,.5477/.4461], [ 15, .5094/.4461] ], 'cptb': [[ -15,.4856/.4461], [ 15,.4839/.4461] ], 'ctW': [[ -2, .5820/.4461], [ 2, .6182/.4461] ], 'ctZ': [[ -2,.5617/.4461], [ 2,.5439/.4461] ],
        'ctp': [[ -11,1.264/.4461], [ 40,1.016/.4461] ]
    },
    'ttbb' : {
        'SM': [[0,1]], 
        'cbW':[[-7,9.975/8.763],[7,10.0/8.763]], 'cpQ3': [[ -7, 8.375/8.763], [ 7, 10.64/8.763] ], 'cpQM': [[ -12, 8.763/8.763], [ 16,8.763/8.763] ],
        'cpt': [[ -20,8.762/8.763], [ 15, 8.764/8.763] ], 'cptb': [[ -15, 9.178/8.763], [ 15,9.262/8.763] ], 'ctW': [[ -2, 8.914/8.763], [ 2, 8.853/8.763] ], 'ctZ': [[ -2, 8.851/8.763], [ 2, 8.774/8.763] ],
        'ctp': [[ -11, 8.762/8.763], [ 40, 8.763/8.763] ]
    },
}

wc_ranges = {
    #'ctG'  :np.arange(-1,1,.01) 
    'cbW': np.arange(-14,14,1),
    'cpQ3': np.arange(-14,14,1),
    'cpQM': np.arange(-20,24,1),
    'cpt': np.arange(-40,40,2),
    'cptb': np.arange(-30,30,2),
    'ctW': np.arange(-4,4,.2),
    'ctZ': np.arange(-4,4,.2),
    'ctp': np.arange(-30,60,3),
}

def main(input_files, samples):
    norm_change = {}
    for input_file, sample in zip(input_files, samples):
        print (input_file, sample)
        norm_change[sample] = process(input_file,re.sub(r'_201\d','',sample))
        #plot_overal_norm(norm_change, sample)
    #
    save_pdf(re.sub(r'_201\d','',samples[0])+"_incxs_impacts.pdf")(plot_overall_norm)(norm_change, samples)
    #plot_overall_norm_compare(norm_change, samples)


def process(input_file, sample):
    eft_df = pd.DataFrame()
    with open(input_file) as ifile:
        eft_df = __worker(ifile.readlines())

    # get parameterization
    df = getBeta(sample).calcBeta(eft_df,sample)
    return get_norm_change(df, sample)
    #plot_overal_norm(norm_change, sample)
    
def get_norm_change(df,sample):
    norm_change = {}
    #pool = Pool()
    for wc, r in wc_ranges.items():
        _norm = partial(calc_norm, df=df, wc=wc)
        #norm_change[wc] = [calc_norm(df,wc,v) for v in r]
        norm_change[wc] = [_norm(v) for v in r]
        #norm_change[wc] = map(_norm, r)
        #pool.close()
    return norm_change

@save_pdf("ttH_2018_incxs_impacts.pdf")
def plot_overall_norm_compare(norm_change, samples):
    for w in norm_change[samples[0]]: # samples[0] is dummy index
        fig, ax = plt.subplots(1,len(samples), sharey=True) # anticipating this is ttz vs tth
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.93,
            hspace=0.0,
            wspace=0.2
        )
        test_points = { sample: np.array(single_points[sample]['SM']+single_points[sample][w]) for sample in samples}
        for i, sample in enumerate(samples):
            ax[i].plot(wc_ranges[w],norm_change[sample][w])
            ax[i].scatter(x=test_points[sample][:,0], y=test_points[sample][:,1], 
                       marker='*',color='red',label='Validation')
            ax[i].axhline(1,c='r',ls='--')
            #plt.xticks([i for i in range(len(norm_change))], norm_change.keys())
            ax[i].set_xlabel(f'WC {w} ')
            ax[i].set_ylabel(rf'$\sigma$({sample})'+r'$^{EFT}$/$\sigma$'+f'({sample})'+r'$^{SM}$')
            ax[i].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i].yaxis.set_minor_locator(AutoMinorLocator())
            ax[i].set_title(sample)
            ax[i].grid(True)
        #fig.suptitle(f'Inc. Rate change w/resp. to WC {w}}')
        #plt.xlim(-1,1)
        #plt.legend()
        #plt.show()



def plot_overall_norm(norm_change, samples):
    for sample in samples:
        nc = norm_change[sample]
        for w in nc:
            test_points = np.array(single_points[re.sub(r'_201\d','',sample)]['SM']+single_points[re.sub(r'_201\d','',sample)][w])
            fig, ax = plt.subplots()
            ax.plot(wc_ranges[w],nc[w])
            ax.scatter(x=test_points[:,0], y=test_points[:,1], marker='*',color='red',label='Validation')
            ax.axhline(1,c='r',ls='--')
            #plt.xticks([i for i in range(len(norm_change))], norm_change.keys())
            ax.set_ylim(0.5,max(ax.get_ylim()[1],1.5))
            ax.set_xlabel(f'WC {w} ')
            ax.set_ylabel(rf'$\sigma$({sample})'+r'$^{EFT}$/$\sigma$'+f'({sample})'+r'$^{SM}$')
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            fig.suptitle(f'Inc. Rate change w/resp. to {w} for {sample}')
            #plt.xlim(-1,1)
            #plt.legend()
            #plt.show()

def __worker(lines):
    pool = Pool(8)
    results = pool.map(get_eft_df,[line.strip('\n') for line in lines])
    pool.close()
    _out_df = pd.concat(results, axis='rows', ignore_index=True)
    return _out_df

def get_eft_df(roofile):
    print(roofile)
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

def calc_norm(v,df=None,wc=None):
    p = sum(df[f'{wc}_{wc}']*v*v)
    q = sum(df[f'{wc}']*v)
    r = sum(df['SM'])
    return p/r + q/r + r/r


if __name__ == '__main__':
    #main("ttjets_eft.txt", 'ttjets')
    input_files = [
        #"files/ttz_eft_2018.txt",
        
        #"files/tth_eft_2016.txt",
        #"files/tth_eft_2017.txt",
        #"files/tth_eft_2018.txt",
        #"files/ttz_eft_2016.txt",
        #"files/ttz_eft_2017.txt",
        #"files/ttz_eft_2018.txt",
        
        "files/tt_B_eft_2018.txt"
    ]
    #main(input_files, ['ttZ','ttH'])
    #main(input_files, ['ttH_2016','ttH_2017','ttH_2018'])
    #main(input_files, ['ttZ_2017','ttZ_2018'])
    main(input_files, ['ttbb_2018'])
    #main(input_files, ['ttbb'])

