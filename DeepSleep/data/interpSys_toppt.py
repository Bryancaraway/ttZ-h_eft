import uproot
import os
import sys
sys.path.insert(1,'/home/bcaraway/ttZh_ana/DeepSleep')
from modules.AnaDict import AnaDict
import numpy as np


def main():
    roofile = 'LostLepton_topPt_systematics.root'
    h_list = ['topPt_up','topPt_dn']
    
    sf = AnaDict()
    if os.path.exists(roofile):
        with uproot.open(roofile) as f_:
            for h in h_list:
                if h in f_:
                    print(roofile, h)
                    sf[h] = {}
                    hist = f_[h]
                    pt_bins = hist.edges
                    v    = np.array(hist.values)
                    verr = hist.variances
                    for i in range(v.shape[0]):
                        sf[h][f'{pt_bins[i]},{pt_bins[i+1]}'] = v[i]
                    #
                #
            #
        #
        sfdir    = '/home/bcaraway/ttZh_ana/DeepSleep/files/toppt_sys_files/'
        sf_out_file = 'toppt_sys.pkl'
        sf.to_pickle(sfdir+sf_out_file)

if __name__ == '__main__':
    #main('electron')
    main()
                        
