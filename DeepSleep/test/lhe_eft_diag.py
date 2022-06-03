import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import uproot
from uproot_methods import TLorentzVectorArray 
getTLVm = TLorentzVectorArray.from_ptetaphim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from functools import partial
#
import config.ana_cff as cfg
from modules.AnaDict import AnaDict
from eftnorm_validation import getBeta
from lib.fun_library import save_pdf

pdgId = {
     25:'H' ,
     23:'Z' ,
     1: 'd' ,
     2: 'u' ,
     3: 's' ,
     4: 'c' ,
     5: 'b' ,
     6: 't' ,
     21:'g' ,
     0: '0',
    11: 'el',
    12: 'enu',
    13: 'mu',
    14: 'mnu',
    15: 'tau',
    16: 'tnu',
}
sig_pdgId = {'H':25, 'Z':23}

part_prod = {'Z':'ttZ','H':'ttH','b':'ttbb'}
wcs = {'ctZ':1.11,'ctW':0.96,'cpQ3':3.04,'cpQM':8.73,'cpt':6.33,'ctp':30.10,'cptb':10.75,'cbW':4.18}

def sort_eft_contrib_by_eft(f_names, f_part):
    wfunc = partial(__worker, f_part=f_part)
    results = map(wfunc,f_names)
    eft_df = pd.concat(results, axis='rows', ignore_index=True)
    #
    for pdf_extra in ['pdf','extra']:
        for wc in wcs:
            fig, ax = plt.subplots()
            fig.suptitle(f'{pdf_extra} parts, {f_part} pT [300,450] GeV')
            ind = np.arange(len(np.unique(eft_df[f'{pdf_extra}_parts'])))
            #eft_norms = [calc_norm(wcs[wc],eft_df[eft_df[f'{pdf_extra}_parts'] == part],wc) for part in np.unique(eft_df[f'{pdf_extra}_parts'])] 
            eft_norms = np.array([calc_yield(wcs[wc],eft_df[eft_df[f'{pdf_extra}_parts'] == part],wc) for part in np.unique(eft_df[f'{pdf_extra}_parts'])])
            SM_norms =  np.array([calc_yield(0,eft_df[eft_df[f'{pdf_extra}_parts'] != part],wc) for part in np.unique(eft_df[f'{pdf_extra}_parts'])])
            eft_norms = (eft_norms + SM_norms) / (calc_yield(0, eft_df, wc))
            eft_yields = np.array([calc_yield(wcs[wc],eft_df[eft_df[f'{pdf_extra}_parts'] == part],wc) for part in np.unique(eft_df[f'{pdf_extra}_parts'])])/(calc_yield(wcs[wc], eft_df, wc))
            ax.scatter(ind,eft_norms,label=part_prod[f_part]+' '+wc+f' = {wcs[wc]}')
            ax.axhline(calc_norm(wcs[wc],eft_df,wc), c='r', label=f'total impact ({len(eft_df)})')
            ax.axhline(1, c='k', label='SM')
            ax.set_xticks(ind)
            ax.set_xticklabels(zip(np.unique(eft_df[f'{pdf_extra}_parts']),[f'{100*y:.2f}%' for y in eft_yields]), rotation =45)
            #ax.set_ylabel('(EFT+SM)/SM')
            ax.set_ylabel('SMEFT(part)+SM(others) / SM(all)')
            #ax.legend(edgecolor='none')
            ax.legend()
            plt.tight_layout()
            #plt.show()
            #exit()

def __worker(f_name, f_part):
    tree = uproot.open(f_name)['Events']
    lhe_parts = AnaDict({ k_: tree.array(k_).pad(6 if f_part == 'H' else 7).fillna(0) for k_ in ['LHEPart_pdgId','LHEPart_pt','LHEPart_eta','LHEPart_phi','LHEPart_mass','LHEPart_status'] })
    eft_w = tree.array('LHEReweightingWeight').pad(184).fillna(1)
    eft_weights = pd.DataFrame.from_dict({f'EFT{i}': eft_w[:,i] for i in range(184)})
    eft_df = getBeta(part_prod[f_part]).calcBeta(eft_weights,part_prod[f_part]).filter(regex=r'c|SM')
    # 
    boosted_cut = get_boosted_cut(lhe_parts,f_part).flatten()
    lhe_parts = lhe_parts[boosted_cut]
    eft_df  = eft_df[boosted_cut]
    # add initiating parts to df
    eft_df['pdf_parts'] = [''.join(sorted([pdgId[abs(row[0])],pdgId[abs(row[1])]]))  for row in lhe_parts['LHEPart_pdgId'][lhe_parts['LHEPart_status'] == -1]]
    eft_df['extra_parts'] = [pdgId[abs(row[-1])] for row in lhe_parts['LHEPart_pdgId']]
    return eft_df
    
def get_boosted_cut(p_, f_):
    if f_ == 'H':
        sig_part = p_['LHEPart_pt'][abs(p_['LHEPart_pdgId']) == sig_pdgId[f_]]
        return ((sig_part > 300) & (sig_part < 450))
    elif f_ == 'Z': # size 7, 5-6 are Z products
        z1, z2 = [getTLVm(p_['LHEPart_pt'][:,i_], p_['LHEPart_eta'][:,i_], p_['LHEPart_phi'][:,i_], p_['LHEPart_mass'][:,i_]) for i_ in [5,6]]
        return (((z1+z2).pt > 300) & ((z1+z2).pt < 450))
        
def calc_norm(v,df=None,wc=None):
    p = sum(df[f'{wc}_{wc}']*v*v)
    q = sum(df[f'{wc}']*v)
    r = sum(df['SM'])
    return p/r + q/r + r/r

def calc_yield(v,df=None,wc=None):
    return calc_norm(v,df,wc)*sum(df['SM'])

def calc_norm_trunc(df=None,wc=None):
    p = df[f'{wc}_{wc}']
    q = df[f'{wc}']
    r = df['SM']
    out = p/r + q/r + r/r
    return out[(out<np.percentile(out, 97.5)) & (out>np.percentile(out, 2.5))]

if __name__ == '__main__':
    f_header = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2016/'
    f_names = [f_header+ 'TTH_EFT_2016/' + f_ for f_ in ''' 
    prod2016MC_v7_NANO_110_Skim_79.root   prod2016MC_v7_NANO_112_Skim_168.root  prod2016MC_v7_NANO_114_Skim_25.root   prod2016MC_v7_NANO_116_Skim_171.root  
    prod2016MC_v7_NANO_118_Skim_34.root   prod2016MC_v7_NANO_11_Skim_136.root prod2016MC_v7_NANO_111_Skim_22.root   prod2016MC_v7_NANO_113_Skim_198.root  
    prod2016MC_v7_NANO_115_Skim_5.root    prod2016MC_v7_NANO_117_Skim_37.root   prod2016MC_v7_NANO_119_Skim_2.root '''.split()]
    save_pdf('lhe_ttH_eft_impact_parts.pdf')(sort_eft_contrib_by_eft)(f_names, 'H')
    f_names = [f_header+ 'TTZ_EFT_2016/' + f_ for f_ in ''' 
    prod2016MC_v7_NANO_110_Skim_182.root  prod2016MC_v7_NANO_112_Skim_192.root  prod2016MC_v7_NANO_114_Skim_94.root   prod2016MC_v7_NANO_116_Skim_159.root  
    prod2016MC_v7_NANO_118_Skim_39.root   prod2016MC_v7_NANO_11_Skim_199.root   prod2016MC_v7_NANO_111_Skim_64.root   prod2016MC_v7_NANO_113_Skim_76.root   
    prod2016MC_v7_NANO_115_Skim_89.root   prod2016MC_v7_NANO_117_Skim_187.root  prod2016MC_v7_NANO_119_Skim_100.root '''.split()]
    save_pdf('lhe_ttZ_eft_impact_parts.pdf')(sort_eft_contrib_by_eft)(f_names, 'Z')
