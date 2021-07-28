import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
import numpy as np
import re
import seaborn as sns
import config.ana_cff as cfg
from lib.fun_library import save_pdf, getLaLabel
from qcDatacard import tbins_map
from post_fit import PostFit
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib import rc

rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            

fit_vars = cfg.withbbvl_dnn_ZHgenm_vars
#fit_vars = [f'NN_{i}' for i in range(64)]

@save_pdf('qc_nn_postfits.pdf')
#@save_pdf('qc_nn_postfits_hl2.pdf')
#@save_pdf('qc_zhm_postfit.pdf')
def main():
    print(len(fit_vars))
    for v in fit_vars:
        print(v)
        froo = f'fitDiagnostics_{v}_NNcuts_run2.root'
        qc= QCNNPostFit(froo,v,tbins_map[v])
        qc.makeplots(doPull=True)

class QCNNPostFit(PostFit):
    '''
    steal some of the functionality 
    from already existing class
    '''

    def __init__(self, fitroo, kinem, t_labels):
        super().__init__(fitroo, kinem=kinem, t_labels=t_labels)
        
    def do_stack(self, d, top_axs, ch):

        ycumm = None
        # stack portion
        #ordered_list = re.findall(rf'tt[H,Z]', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_2b','tt_bb','TTBar']
        ordered_list = re.findall(rf'tt[H,Z]', ' '.join(list(d.keys()))) + ['single_t','VJets','ttX','tt_B','TTBar']
        #colors =  plt.cm.gist_rainbow(np.linspace(0,1,len(ordered_list)))
        colors =  plt.cm.tab10(np.linspace(0,1,10))[0:2]
        colors = np.append(colors, plt.cm.gist_rainbow(np.linspace(0,1,6)), axis=0)
        #
        #for j,k in enumerate(d):
        for j,k in enumerate(ordered_list):
            #if 'data' in k or 'total' in k: continue
            if k not in d: continue
            y = np.append(d[k]['values'],0)
            if ycumm is not None:
                ycumm += y
            else:
                ycumm = y 
            c = colors[j]
            #c = colors[j + (len(colors)//2)*(j % 2) - 1*((j+1)%2)]
            label = getLaLabel(k)[0]
            top_axs.fill_between(self.edges[ch],ycumm,ycumm-y,step="post", 
                                 linewidth=0, color=c, label=label)
            # add total error and data points
        top_axs.errorbar(x=(self.edges[ch][1:]+self.edges[ch][:-1])/2, y=d['data']['values'],
                     xerr=(self.edges[ch][1:]-self.edges[ch][:-1])/2 ,yerr=[d['data']['errdw'],d['data']['errup']], 
                     fmt='.', label='data', color='k')
        self.make_error_boxes(top_axs, (self.edges[ch][1:]+self.edges[ch][:-1])/2, d['total']['values'],
                              xerror=(self.edges[ch][1:]-self.edges[ch][:-1])/2, yerror=d['total']['err'], label='stat+sys')
    


if __name__ == '__main__':
    main()
