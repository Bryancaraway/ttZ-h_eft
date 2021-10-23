import sys
import os
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
import numpy as np


#DIR = "/cms/data/hatake/ana/TTX/Powheg-CMS/CMSSW_10_6_20/src/genproductions/bin/Powheg/"
#SUB_DIRS = ["run_ttH_default","run_ttH_fakevirt0","run_ttH_fixedscale"]
DIR = "/cms/data/hatake/ana/TTX/Powheg-CMS/CMSSW_9_3_0/src/genproductions/bin/Powheg/"
SUB_DIRS = ["run_ttH_v2"]
LHE_FILE = "cmsgrid_final.lhe"


def main():
    fig, ax = plt.subplots()
    # open files and store info
    for sub_dir in SUB_DIRS:
        with open(f"{DIR}{sub_dir}/{LHE_FILE}") as lhe_file:
            lines = lhe_file.readlines()
            get_scale_weights(lines, sub_dir.lstrip("run_ttH_"), ax)
    #plt.xlim(0,2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.legend(ncol=2)
    plt.show()

def get_scale_weights(lines, ftype, ax):
    w_strs_nom = re.findall(r"<wgt id='1001'>[-]?\d*[.,]?\d*E[+-]?\d*</wgt>", " ".join(lines))
    w_strs_up  = re.findall(r"<wgt id='1005'>[-]?\d*[.,]?\d*E[+-]?\d*</wgt>", " ".join(lines))
    w_strs_dn  = re.findall(r"<wgt id='1009'>[-]?\d*[.,]?\d*E[+-]?\d*</wgt>", " ".join(lines))
    #
    w_nom = np.array([float(w_str) for w_str in re.findall(r"[-]?\d*[.,]?\d*E[+-]?\d*" , " ".join(w_strs_nom))])
    w_up  = np.array([float(w_str) for w_str in re.findall(r"[-]?\d*[.,]?\d*E[+-]?\d*" , " ".join(w_strs_up))])
    w_dn  = np.array([float(w_str) for w_str in re.findall(r"[-]?\d*[.,]?\d*E[+-]?\d*" , " ".join(w_strs_dn))])

    print(len(w_nom),len(w_up),len(w_dn))
    ax.hist(w_up/w_nom, bins=20, range=(.5,1.5), histtype="step", label=f"{ftype} Scale Up")
    ax.hist(w_dn/w_nom, bins=20, range=(.5,1.5), histtype="step", label=f"{ftype} Scale Down")

    

    

if __name__ == '__main__':
    main()
