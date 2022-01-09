import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import numpy as np
import pandas as pd
import json
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
from modules.AnaDict import AnaDict
from lib.fun_library import save_pdf, clop_pear_ci, getFakebbvlCuts
import config.ana_cff as cfg

impacts_dir  = 'fitdiag_roots/'
impacts_file = 'impacts_run2.incv2.json'

top_nuisances = 8
sig = {'r_ttH' : 0.86, 'r_ttZ': 1.02}

def calc_pdiff():
    i_file = json.load(open(impacts_dir+impacts_file,'r'))
    i_dict = {param_dict['name'] : {'r_ttH':param_dict['impact_r_ttH'], 'r_ttZ':param_dict['impact_r_ttZ']} for param_dict in i_file['params'] if 'prop' not in param_dict['name']}
    top_tth = dict(sorted(i_dict.items(), key=lambda item: -item[1]['r_ttH']))
    top_ttz = dict(sorted(i_dict.items(), key=lambda item: -item[1]['r_ttZ']))
    nuisances = set([list(top_tth)[i] for i in range(top_nuisances)] + [list(top_ttz)[i] for i in range(top_nuisances)])
    print(f"{'var':17} | {'ttZ':3} | {'ttH':3} ")
    for nui in nuisances:
        print(f"nui: {nui:12} | {pdiff_form(i_dict[nui], 'r_ttZ'):.1f} | {pdiff_form(i_dict[nui], 'r_ttH'):.1f}")


def pdiff_form(_x, _k):
    dsig = (np.abs(sig[_k] - np.sqrt( np.power(sig[_k],2) - np.power(_x[_k],2))) / sig[_k]) * 100
    return dsig

if __name__ == '__main__':
    calc_pdiff()
