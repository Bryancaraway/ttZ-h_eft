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
import uproot
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import config.ana_cff as cfg




class PostFit:
    '''
    script that makes post-fit and pre-fit 
    stack plots of data vs MC from fitDiagnostic file
    '''
    #
    hists = {} # dictionary with format Xfit:ZhptX:process:values
    edges = {} # dictionary with format ZhptX:edges
    
    def __init__(self,fitroo):
        self.fitroo = fitroo
        self.year = re.search(r'(201\d|run2)',fitroo).group()
        self.load_info()

    def load_info(self):
        with uproot.open(self.fitroo) as roo:
            # store, prefit # shapes_prefit, shapes_fit_s
            def to_dict(p): # takes prefit or postfit
                pp_dict = {'prefit' :'shapes_prefit',
                           'postfit':'shapes_fit_s'}
                suf = (lambda b: b.decode().split(';')[0])
                self.hists[p] = {}
                for pt in roo[pp_dict[p]]:
                    pt_str = suf(pt)
                    self.hists[p][pt_str] = {}
                    if pt_str not in self.edges:
                        self.edges[pt_str]=roo[pp_dict[p]][pt]['TTBar'].edges
                    for hist in roo[pp_dict[p]][pt]:
                        hist_str = suf(hist)
                        h = roo[pp_dict[p]][pt][hist]
                        try:
                            self.hists[p][pt_str][hist_str] = {'values':h.values,  'err':np.sqrt(h.variances)}
                        except:
                            self.hists[p][pt_str][hist_str] = {'values':h.yvalues, 'errup':h.yerrorshigh,'errdw':h.yerrorslow}
            #
            to_dict('prefit')
            to_dict('postfit')
        #
    #
    def makeplots(self):
        self.init_axes()
        self.make_stackhist()
    #
    def init_axes(self):
        # init fig and four axes
        #fig, (ax1t, ax1b, ax2t, ax2b) = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig, axs = plt.subplots(2,2, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.0,
            wspace=0.2
        )
        self.fig = fig
        self.axs = np.array(axs) # [1:[t,b],2:[t,b]]
    #
    def make_stackhist(self):
        ax = self.axs[0,0]
        test = self.hists['prefit']['Zhpt1']
        for k in test:
            plst.step()
        #plt.stackplot(*[test[k]['values'] for k in test if k != 'data' or 'total' not in k])
        #plt.show()

if __name__ == '__main__':
    # fitDiagnostics_test_2016.root, fitDiagnostics_test_2017.root, fitDiagnostics_test_2018.root, sfitDiagnostics_test_run2.root
    PostFit("fitDiagnostics_test_2016.root").makeplots()
