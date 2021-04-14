import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import matplotlib.pyplot as plt
from matplotlib import rc
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

import json
import numpy as np
import config.ana_cff as cfg
from lib.fun_library import getZhbbBaseCuts, getZhbbWeight, save_pdf

mfp = cfg.master_file_path
ana_cuts = getZhbbBaseCuts

def make_2deffsfplot(lep,year): # "Electron/Muon"
    eff_dict = json.load(open(cfg.dataDir+'/lep_effsf_files/'+f"trigeffSF_{lep}_{year}.json", 'r'))[lep]
    fig, ax = plt.subplots()
    fig.subplots_adjust(
        top=0.88,
        bottom=0.11,
        left=0.11,
        right=0.99,
        #hspace=0.0 
        wspace=0.0
        )
    c = ax.pcolor(eff_dict['pt_bins'],eff_dict['eta_bins'],np.array(eff_dict['pt_eta_sf']).T, vmin=0, vmax=1.05)
    #c = ax.pcolor(eff_dict['pt_eta_sf'])
    for i in range(len(eff_dict['eta_bins'])-1):
        for j in range(len(eff_dict['pt_bins'])-1):
            ax.text( 
                eff_dict['pt_bins'][j] + (eff_dict['pt_bins'][j+1]-eff_dict['pt_bins'][j])/2, 
                eff_dict['eta_bins'][i] + (eff_dict['eta_bins'][i+1]-eff_dict['eta_bins'][i])/2, 
                f"{eff_dict['pt_eta_sf'][j][i]:.3f}\n"+r"${{}}^{{+{0:.3f} ({1:.3f},{4:.3f})}}_{{-{2:.3f} ({3:.3f},{4:.3f})}}$".format(
                    eff_dict['pt_eta_sf_Up'][j][i],eff_dict['pt_eta_sf_stat_Up'][j][i],
                    eff_dict['pt_eta_sf_Down'][j][i],eff_dict['pt_eta_sf_stat_Down'][j][i],
                    eff_dict['pt_eta_sf_sys'][j][i]),
                horizontalalignment='center', verticalalignment='center', fontsize=4.0)
    #
    ax.set_xscale("Log")
    ax.set_xticks(eff_dict['pt_bins'])
    ax.set_xticklabels([str(i) for i in eff_dict['pt_bins']])
    ax.set_yticks(eff_dict['eta_bins'])
    ax.set_yticklabels([f"{i:.1f}" for i in eff_dict['eta_bins']])
    plt.minorticks_off()
    fig.suptitle(f'{lep}_{year}')
    cbar = fig.colorbar(c, pad=0)
    cbar.set_ticks(np.arange(0,1.1,.1))
    cbar.set_ticklabels([f"{i:.1f}" for i in np.arange(0,1.1,.1)])
    

@save_pdf("leptrig_effsf.pdf")
def main():
    for lep in ['Electron','Muon']:
        for year in cfg.Years:
            make_2deffsfplot(lep,year)    

if __name__ == '__main__':
    main()
