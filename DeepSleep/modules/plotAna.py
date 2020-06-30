###### Plotter code for analysis ######
import pandas as pd
import numpy as np
import re
from functools import partial, reduce
import operator as op
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
rc("savefig",dpi=250)
rc("figure", figsize=(6, 6*(6./8.)), dpi=200)                                                            
#rc("text.latex", preamble=r"\usepackage{amsmath}")                                                             
#rc('font',**{'family':'serif','serif':['Times']})
#rc("hatch", linewidth=0.0) 
#
from lib.fun_library import getZhbbBaseCuts, getLaLabel
import cfg.deepsleepcfg as cfg

class Plotter :
    '''
    master class designed to contain other plotting classes
    main functions: grab files by which to analyze
    setup plotting dpendencies
    '''
    cut_func = staticmethod(getZhbbBaseCuts) # M 50-200, pt 200, b_out_Zh 2
    fontsize = 12
    saveDir  = 'pdf/'
    pklDir   = 'files/'
    data_samples = ['MuData','EleData']
    year      = None
    data_dict = None
    w_dict    = None
    i_dict    = None
    real_data = None
    HEM_opt   = ''

    def __init__(self,samples, kinem, bin_range, 
                 xlabel=None, n_bins=20, 
                 doLog=True, doNorm=False, 
                 doShow = True, doSave = False,
                 doCuts=True,add_cuts=None,sepGenOpt=None,addData=False):

        self.samples    = samples
        self.kinem      = kinem
        self.bin_range  = bin_range
        self.xlabel     = kinem if xlabel is None else xlabel
        self.n_bins     = n_bins
        self.bins       = np.arange(bin_range[0],bin_range[-1]+((bin_range[-1]-bin_range[0])/n_bins) , (bin_range[-1]-bin_range[0])/n_bins)
        self.doLog      = doLog
        self.doNorm     = doNorm
        self.doShow     = doShow
        self.doSave     = doSave
        self.addData    = addData
        self.doCuts     = doCuts
        self.add_cuts   = add_cuts if add_cuts is not None else ''
        self.sepGenOpt  = sepGenOpt
        #
        self.prepData()

    def prepData(self):
        data = {sample: self.apply_cuts(self.data_dict[sample]) for sample in self.samples}
        if self.addData:
            data_dir = f'{self.pklDir}{self.year}/data_files/'
            self.real_data  = {'Data': np.clip(np.hstack(
                [self.apply_data_cuts(sample,pd.read_pickle(f'{data_dir}{sample}_val.pkl'))[self.kinem].values 
                 for sample in self.data_samples]),self.bin_range[0],self.bin_range[-1])}
        #
        self.data = data
        self.sepData()
        print(self.year+self.HEM_opt)
        #
        self.w_dict = {k: v['weight']* np.sign(v['genWeight']) 
                       * (self.lumi/cfg.Lumi[self.year])
                       * v['Stop0l_topptWeight']
                       * (v['SAT_HEMVetoWeight_drLeptonCleaned']  if self.year+self.HEM_opt == '2018' else 1.0 )
                       #* (v['Stop0l_topMGPowWeight'] if self.year == '2017' else 1.0)
                       * (pd.concat([v['Stop0l_trigger_eff_Electron_pt'][v['passSingleLepElec']==1],v['Stop0l_trigger_eff_Muon_pt'][v['passSingleLepMu']==1]]).sort_index() if self.year != '2015' else 1.0)
                       * v['BTagWeight'] 
                       * v['puWeight']  
                       * (v['PrefireWeight'] if self.year != '2018' else 1.0)
                       * v['ISRWeight']  
                       for k,v in self.data.items()}

        self.i_dict = {k: sum(v) for k,v in self.w_dict.items()}
        #
        self.data = {k: np.clip(v[self.kinem],self.bin_range[0],self.bin_range[-1]) for k,v in self.data.items()}
    
    def apply_data_cuts(self, leptype, df):
        temp_add_cuts = self.add_cuts
        add_by_lep = {'EleData':'passSingleLepElec==1',
                      'MuData' :'passSingleLepMu==1'}
        add_by_HEM = {'preHEM' :';run<319077',
                      'postHEM':';run>=319077',
                      '':''}
        data_cuts = add_by_lep[leptype]+add_by_HEM[self.HEM_opt]
        if self.add_cuts is '':
            self.add_cuts = data_cuts
        else:
            self.add_cuts += ';'+data_cuts
        df = self.apply_cuts(df)
        self.add_cuts = temp_add_cuts
        return df

    def apply_cuts(self,df):
        df_mask = self.cut_func(df) if self.doCuts else True
        if self.add_cuts is '': 
            if self.doCuts:
                return df[df_mask]
            else:
                return df

        b_str    = re.findall(r'[><=!]=*',self.add_cuts)
        interp_dict = {'>=':op.ge, '<=':op.le, '>':op.gt, '<':op.lt, '==':op.eq, '!=':op.ne}
        for cut, k in zip(self.add_cuts.split(';'),b_str):
            kinem, val = cut.split(k)
            df_mask = interp_dict[k](df[kinem],float(val)) & df_mask
        return df[df_mask]

    def sepData(self):
        if self.sepGenOpt is None: return
        interp_dict = {'sepGenSig': self.sepGenSig,
                       'sepGenBkg': self.sepGenBkg,
                       'sepGenMatchedSig':self.sepGenMatchedSig,
                       'sepGenMatchedBkg':self.sepGenMatchedBkg}
        
        for k in self.sepGenOpt.split(';'):
            interp_dict[k]()
        
    def sepGenSig(self):
        df = self.data
        df['TTZH_Zbb'] = (lambda x : x[x['Zbb'] == True])(df['TTZH'])
        df['TTZH_Hbb'] = (lambda x : x[x['Hbb'] == True])(df['TTZH'])
        df['TTZH_Zqq'] = (lambda x : x[x['Zqq'] == True])(df['TTZH'])
        self.data.update(df)
        self.data.pop('TTZH')

    def sepGenBkg(self):
        df = self.data
        df['TTBarLep_bb']   = (lambda x : x[x['tt_bb'] == True])(df['TTBarLep'])
        df['TTBarLep_nobb'] = (lambda x : x[x['tt_bb'] == False])(df['TTBarLep'])
        self.data.update(df) 
        self.data.pop('TTBarLep')
        
    def sepGenMatchedSig(self):
        df = self.data 
        df['TTZH_genZbb'] = (lambda x : x[x['matchedGen_Zbb'] == True])(df['TTZH'])
        df['TTZH_genHbb'] = (lambda x : x[x['matchedGen_Hbb'] == True])(df['TTZH'])
        df['TTZH_genZqq'] = (lambda x : x[x['matchedGen_Zqq'] == True])(df['TTZH'])
        df['TTZH_noGenMatch'] = (lambda x : x[(x['matchedGen_Zqq'] == False) & 
                                                          (x['matchedGen_ZHbb'] == False)])(df['TTZH'])
        self.data.update(df) 
        self.data.pop('TTZH')

    def sepGenMatchedBkg(self):
        df = self.data
        df['TTBarLep_bbGenMatch']   = (lambda x : x[x['genMatched_tt_bb'] == True])(df['TTBarLep'])
        df['TTBarLep_nobbGenMatch'] = (lambda x : x[x['genMatched_tt_bb'] == False])(df['TTBarLep'])
        self.data.update(df) 
        self.data.pop('TTBarLep')

    @property
    def retData(self):
        k_    = np.array([k for k in self.data])
        h_   = np.array([self.data[k].to_numpy()    for k in self.data])
        w_   = np.array([self.w_dict[k].to_numpy()  for k in self.data])
        i_   = np.array([self.i_dict[k]             for k in self.data])
        l_,c_ = np.array([np.array(getLaLabel(k))   for k in self.data]).T
        l_   = np.array([f'{x} ({y:3.1f})' for x,y in zip(l_,i_)])
        #print( h_,w_,i_,c_,l_)
        return k_,h_,w_,i_,c_,l_

    @property
    def beginPlt(self):
        if self.addData:
            self.fig, (self.ax, self.ax2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
            #self.fig.subplots_adjust(hspace=0)

        else:
            self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.0 if self.addData else 0.2,
            wspace=0.2
        )

    @property
    def endPlt(self):
        self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.tick_params(which='both', direction='in', top=True, right=True)
        self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Simulation$", fontsize = self.fontsize)
        self.fig.text(0.635,0.89, f'{self.lumi}'+r' fb$^{-1}$ (13 TeV)',  fontsize = self.fontsize)
        plt.xlabel(self.xlabel, fontsize = self.fontsize)
        self.ax.set_ylabel(f"{'%' if self.doNorm else 'Events'} / {(self.bin_range[-1]-self.bin_range[0])/self.n_bins:.2f}")#fontsize = self.fontsize)
        plt.xlim(self.bin_range)
        if self.doLog: self.ax.set_yscale('log')
        plt.grid(True)
        #plt.setp(patches_, linewidth=0)
        self.ax.legend(framealpha = 0, ncol=2, fontsize='xx-small')
        if self.doSave: plt.savefig(f'{self.saveDir}{self.xlabel}_.pdf', dpi = 300)
        if self.doShow: plt.show()
        plt.close(self.fig)

    @classmethod
    def load_data(cls,year='2017',HEM_opt=''):
        cls.year = year
        cls.HEM_opt = HEM_opt
        cls.lumi = cfg.Lumi[year+HEM_opt]
        cls.data_dict = {sample: pd.read_pickle(f'{cls.pklDir}{year}/mc_files/{sample}_val.pkl') 
                         for sample in (cfg.MC_samples if year != '2018' else cfg.MC_samples)}
        

class StackedHist(Plotter) :
    #Creates stacked histograms
    
    def __init__(self,samples, kinem, bin_range, 
                 xlabel=None, n_bins=20, 
                 doLog=True, doNorm=False, 
                 doShow = True, doSave = False,
                 doCuts=True,add_cuts=None,sepGenOpt=None,addData=False):
        super().__init__(samples,kinem,bin_range,xlabel,n_bins,doLog,doNorm, 
                         doShow,doSave,doCuts,add_cuts,sepGenOpt,addData)
        #
        self.makePlot()

    def makePlot(self):
        self.beginPlt
        #
        h, w, i, c, l = self.sortInputs()
        n_, bins_, patches_ = self.ax.hist(
            h,
            bins=self.bins, 
            stacked=True,# fill=True,
            #range=range_,
            histtype='stepfilled',
            #density=False,
            #linewidth=0,
            weights = w if not self.doNorm else np.divide(w,i),
            color   = c,
            label   = l
        )
        #
        if self.addData: 
            n_mc, _ = np.histogram(np.hstack(h),bins=self.bins, weights=np.hstack(w),
                                   range=(self.bin_range[0],self.bin_range[1]) )
            self.addDataPlot(n_mc)
        self.endPlt

    def addDataPlot(self,n_mc):
        h = self.real_data['Data']
        n_data,edges = np.histogram(h,bins=self.bins, range=(self.bin_range[0],self.bin_range[1]))
        bin_c = (edges[1:]+edges[:-1])/2
        bin_w = edges[1]-edges[0]
        #
        n_data = n_data #* (137/41.9) # to scale to full run2
        #
        self.ax.errorbar(x=bin_c, y = n_data if not self.doNorm else np.divide(n_data,np.sum(n_data)), 
                         xerr=bin_w/2, yerr=np.sqrt(n_data),
                         fmt='.',  color='k',
                         label=f'Data ({sum(n_data)})')
        #
        y = n_data/n_mc
        yerr = y*np.sqrt(np.power(np.sqrt(n_data)/n_data,2)+np.power(0/n_mc,2))
        self.ax2.errorbar(x=bin_c, y = n_data/n_mc,
                          xerr=bin_w/2, yerr=yerr,
                          fmt='.', color='k')
        self.ax2.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        self.ax2.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax2.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        self.ax2.yaxis.set_major_locator(FixedLocator([.75,1,1.25]))
        self.ax2.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax2.tick_params(which='both', direction='in', top=True, right=True)
        self.ax2.set_ylim(0.5,1.5)
        self.ax2.set_ylabel('data/MC')

    def sortInputs(self):
        k_,h_,w_,i_,c_,l_ = self.retData
        # insert signal in first slot for stack plot
        sig = np.argwhere(k_=='TTZH')
        sigh, sigw, sigi, sigc, sigl = map(lambda x: x[sig], [h_,w_,i_,c_,l_])
        h_,w_,i_,c_,l_ = map(lambda x: np.delete(x,sig), [h_,w_,i_,c_,l_])
        #
        ind_ = np.argsort(i_)
        h_,w_,i_,c_,l_ = [np.insert(x[ind_],0,y) 
                          for x,y in zip([h_,w_,i_,c_,l_],[sigh, sigw, sigi, sigc, sigl])]
        return h_,w_,i_,c_,l_

    
class Hist (Plotter) :
    # creats hists for 1-to-1 comparison
    def __init__(self,samples, kinem, bin_range, 
                 droptt=False, dropNoGenM=True, dropZqq=False,
                 xlabel=None, n_bins=20, 
                 doLog=False, doNorm=True, 
                 doShow = True, doSave = False,
                 doCuts=True,add_cuts=None,sepGenOpt=None,addData=False):

        super().__init__(samples,kinem,bin_range,xlabel,n_bins,doLog,doNorm, 
                         doShow,doSave,doCuts,add_cuts,sepGenOpt,addData)

        if dropNoGenM: [self.data.pop(k) for k in re.findall(r'\w*_no\w*GenMatch',' '.join(self.data)) if k in self.data]
        if droptt    : self.data.pop('TTBarLep_nobb') if 'TTBarLep_nobb' in self.data else None
        if dropZqq   : [self.data.pop(k) for k in re.findall(r'TTZH_\w*Zqq',      ' '.join(self.data)) if k in self.data]
        #
        self.makePlot()

    def makePlot(self):
        self.beginPlt
        #
        k, h, w, i, c, l = self.retData
        n_, bins_, patches_ = self.ax.hist(
            h,
            bins=self.bins, 
            stacked=False,# fill=True,
            #range=range_,
            histtype='step',
            #density=False,
            #linewidth=0,
            weights = w if not self.doNorm else np.divide(w,i),
            color   = c,
            label   = l
        )
        #
        self.endPlt
    #
#

