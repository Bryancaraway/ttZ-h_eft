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



inc_dir =  '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2018/TTbb_SemiLeptonic_2018/'
inc_file =  '2C7E7937-7810-6544-A12F-1CE64CB72545_Skim_1.root'
boost_file = cfg.master_file_path+'/2018/mc_files/ttbb_val.pkl'

@save_pdf("inc_vs_boost_ttbb_genbbpt.pdf")
def main():
    inc_ttbb_pt = get_inc_info()
    boost_ttbb_pt = get_boost_info()
    #
    fig, ax = plt.subplots()
    ax.hist(
        np.clip(inc_ttbb_pt, 0, 500),   
        bins=20, histtype='step', 
        weights=np.ones_like(inc_ttbb_pt)*1/len(inc_ttbb_pt), range=(0,500),  label='Inclusive')
    ax.hist(
        np.clip(boost_ttbb_pt, 0, 500), 
        bins=20, histtype='step', 
        weights=np.ones_like(boost_ttbb_pt)*1/len(boost_ttbb_pt), range=(0,500),  label='Passing Selection')
    #
    ax.set_xlim(0,500)
    ax.set_ylabel('fraction of yield / bin')
    ax.set_xlabel('bb gen pt [GeV]')
    ax.legend()
    fig.suptitle("tt+bb: extra gen bb pt ")
    #plt.show()

def get_inc_info():
    g = 'GenPart_'
    gen_vars = ['genTtbarId', g+'pdgId',g+'genPartIdxMother', g+'pt',g+'eta',g+'phi', g+'mass']
    getTLVm = TLorentzVectorArray.from_ptetaphim
    with uproot.open(inc_dir+inc_file) as roo:
        t = roo.get('Events')
        gentt_bb , gen_ids, gen_mom, gen_pt, gen_eta, gen_phi, gen_mass= map(t.array, gen_vars)
        ext_bb = (lambda c: c[(gen_mom <= 0) & (abs(gen_ids) == 5)][(gentt_bb % 100) >= 51])
        #print(gentt_bb , gen_ids, gen_mom, gen_pt, gen_eta, gen_phi, gen_mass)
        b1 = getTLVm(*map((lambda b: b.pad(2)[:,0]), list(map(ext_bb,[gen_pt,gen_eta,gen_phi, gen_mass]))))
        b2 = getTLVm(*map((lambda b: b.pad(2)[:,1]), list(map(ext_bb,[gen_pt,gen_eta,gen_phi, gen_mass]))))
        #return (b1+b2).pt[(gentt_bb % 100) >= 51]
        genbb_pt = (b1+b2).pt
        return genbb_pt[genbb_pt > 0.0]
        
def get_boost_info():
    df = pd.read_pickle(boost_file)
    df = df[(df[cfg.nn] >= 0) & (df['process'] == 'tt_B')] # events passing our selections
    return df['ttbb_genbb_pt'].to_numpy()

if __name__ == '__main__':
    main()
