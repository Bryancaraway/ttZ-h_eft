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
from lib.fun_library import save_pdf, getLaLabel


def make_1d_r_ratio_plot(r_ratio,process,year,ax): 

    edges = list(range(4,4+len(r_ratio)))
    label,c = getLaLabel(process)
    ax.scatter(x=edges, y=r_ratio, c=c, marker='+', label=label)
    

@save_pdf("btag_rratio.pdf")
def main():
    for year in cfg.Years:
        r_ratio = json.load(open(cfg.dataDir+'/btagw_r_ratio/'+f"btagw_r_ratio_{year}.json", 'r'))
        fig, ax = plt.subplots()
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.90,
            #hspace=0.0 
            wspace=0.0
        )
        for p in ['ttZ','ttH','TTBar','ttbb','single_t','ttX','VJets']:
            make_1d_r_ratio_plot(r_ratio[p]['r_ratio'][4:],p,year,ax)    
        ax.legend(ncol=2)
        ax.set_ylim(0.75,1.)
        ax.set_ylabel('r ratio')
        ax.set_xlabel('jet multiplicity (AK4)')
        fig.suptitle(f"B-tag r ratio ({year})")
        ax.grid(color='k', linestyle=':', alpha=0.25)

if __name__ == '__main__':
    main()
