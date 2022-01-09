import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
np.random.seed(1)
import uproot
from uproot_methods import TLorentzVectorArray
import config.ana_cff as cfg
from lib.fun_library import save_pdf, getLaLabel
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc


inc_dir =  '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2018/'
# ttbb
#inc_file =  'TTbb_SemiLeptonic_2018/2C7E7937-7810-6544-A12F-1CE64CB72545_Skim_1.root'
#boost_file = cfg.master_file_path+'/2018/mc_files/ttbb_val.pkl'
# ttbar
inc_file = 'TTToSemiLeptonic_2018/7E1D5911-1000-3546-97BC-51870CB9A867_Skim_121.root'
boost_file = cfg.master_file_path+'/2018/mc_files/TTBar_val.pkl'

b_or_c_quark = 'c'
#@save_pdf("inc_vs_boost_ttbb_genbbpt.pdf")
@save_pdf(f"inc_vs_boost_tt{2*b_or_c_quark}_gen{2*b_or_c_quark}pt.pdf")
def main():
    inc_ttx_pt = get_inc_info(b_or_c_quark)
    quark_boost_info = {'b': get_boost_bb_info, 'c': get_boost_cc_info}
    boost_ttx_pt = quark_boost_info[b_or_c_quark]()
    #
    fig, ax = plt.subplots()
    ax.hist(
        np.clip(inc_ttx_pt, 0, 500),   
        bins=20, histtype='step', 
        weights=np.ones_like(inc_ttx_pt)*1/len(inc_ttx_pt), range=(0,500),  label='Inclusive')
    ax.hist(
        np.clip(boost_ttx_pt, 0, 500), 
        bins=20, histtype='step', 
        weights=np.ones_like(boost_ttx_pt)*1/len(boost_ttx_pt), range=(0,500),  label='Passing Selection')
    #
    ax.set_xlim(0,500)
    ax.set_ylabel('fraction of yield / bin')
    ax.set_xlabel(f'{2*b_or_c_quark} gen pt [GeV]')
    ax.legend()
    fig.suptitle(f"tt+{2*b_or_c_quark}: extra gen {2*b_or_c_quark} pt ")
    #plt.show()

def get_inc_info(bc_quark='b'):
    g = 'GenPart_'
    gen_vars = ['genTtbarId', g+'pdgId',g+'genPartIdxMother', g+'pt',g+'eta',g+'phi', g+'mass']
    getTLVm = TLorentzVectorArray.from_ptetaphim
    quark_cut = {'b': (lambda gttid : (gttid % 100) >= 51), 
                 'c': (lambda gttid : ( (gttid>=41) & (gttid<50) ))}
    quark_id = {'b':5, 'c':4}
    with uproot.open(inc_dir+inc_file) as roo:
        t = roo.get('Events')
        gentt_xx , gen_ids, gen_mom, gen_pt, gen_eta, gen_phi, gen_mass= map(t.array, gen_vars)
        ext_xx = (lambda c: c[(gen_mom <= 0) & (abs(gen_ids) == quark_id[bc_quark])][ quark_cut[bc_quark](gentt_xx) ])
        #print(gentt_bb , gen_ids, gen_mom, gen_pt, gen_eta, gen_phi, gen_mass)
        x1 = getTLVm(*map((lambda _x: _x.pad(2)[:,0]), list(map(ext_xx,[gen_pt,gen_eta,gen_phi, gen_mass]))))
        x2 = getTLVm(*map((lambda _x: _x.pad(2)[:,1]), list(map(ext_xx,[gen_pt,gen_eta,gen_phi, gen_mass]))))
        genxx_pt = (x1+x2).pt
        return genxx_pt[genxx_pt > 0.0]
        
def get_boost_bb_info():
    df = pd.read_pickle(boost_file)
    df = df[(df[cfg.nn] >= 0) & (df['process'] == 'tt_B')] # events passing our selections
    return df['ttbb_genbb_pt'].to_numpy()

def get_boost_cc_info():
    df = pd.read_pickle(boost_file)
    df = df[(df[cfg.nn] >= 0) & (df['tt_C'] == True)] # events passing our selections
    print(df['ttcc_gencc_pt'])
    return df['ttcc_gencc_pt'].to_numpy()

if __name__ == '__main__':
    main()
