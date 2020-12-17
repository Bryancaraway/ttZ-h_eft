import uproot
from uproot_methods import TLorentzVectorArray
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def main(roo, iseft=False):
    t = uproot.open(roo)['Events']
    
    if iseft:
        eftw = t.array('LHEReweightingWeight')[:,183]

    gen_st = t.array('GenPart_statusFlags')[:,4:6] # 4481
    gen_id = t.array('GenPart_pdgId')[:,4:6] # 4481
    kinems = ['pt','eta','phi','mass']
    gen_k = {k: t.array(f'GenPart_{k}')[:,4:6] for k in kinems}
    
    getTLVm = TLorentzVectorArray.from_ptetaphim
    
    # only get b-extra quarks +bb portion
    b_k = {}
    cut = (abs(gen_id) == 5) & (gen_st == 4481)
    
    for k in gen_k:
        
        b_k[k] = gen_k[k][(abs(gen_id) == 5) & (gen_st == 4481)][cut.sum() == 2]

    
    b1_tlv = getTLVm(b_k['pt'][:,0],b_k['eta'][:,0],b_k['phi'][:,0],5)#b_k['mass'][:,0])
    b2_tlv = getTLVm(b_k['pt'][:,1],b_k['eta'][:,1],b_k['phi'][:,1],5)#b_k['mass'][:,1])
    
    bb_tlv = b1_tlv + b2_tlv
    #evcut = ((bb_tlv.pt >= 60) & (bb_tlv.pt < 450) & (bb_tlv.mass >= 20) & (bb_tlv.mass < 200))
    evcut = ((bb_tlv.mass >= 0))
    if iseft:
        eftw = eftw[cut.sum() == 2]

    if iseft:
        plt.hist2d(
            bb_tlv.mass,
            eftw.clip(0,1000),
            bins=[np.linspace(0,50,50),np.logspace(np.log10(10**(-8)),np.log10(1000), 100)],
              norm=colors.LogNorm(),
            )
        plt.colorbar()
        plt.yscale('log')
        plt.show()
        #plt.hist(
        #    eftw.clip(0,1000),
        #    range=(0,1000),
        #    bins=np.logspace(np.log10(10**(-8)),np.log10(1000), 100),
        #    histtype='step',
        #)
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.show()

    plt.hist(
        np.clip(bb_tlv.mass[evcut].tolist(), 0, 20), 
        #b1_tlv.delta_r(b2_tlv)[evcut],
        #range=(0,5), 
        range=(0,20), 
        weights=(eftw[evcut]/sum(eftw[evcut]) if iseft else np.ones_like(b_k['pt'][evcut])/sum(evcut)), 
        #weights=eftw[evcut] if iseft else np.ones_like(b_k['pt'][evcut]), 
        bins=20, histtype='step', 
        #density=True,
        label=roo.replace('.root',''))


main('ttbb_test.root', iseft=True)
main('ttbb_cen.root')


#plt.yscale('log')

plt.legend()
plt.show()


