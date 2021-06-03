#                      #  
##                    ##
########################                               
### TTZ/H, Z/H to bb ###
### build datacard   ###                               
### for HCT          ###                               
########################                               
### written by:      ###                               
### Bryan Caraway    ###                               
########################                               
##                    ##                                 
#                      #

import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import os
#import pickle
#import math
import functools
#from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Pool
import re
#import operator as OP
#
import config.ana_cff as cfg
import config.process_norms as p_norms
from lib.fun_library import weighted_quantile, getZhbbBaseCuts, getZhbbWeight, t2Run
from lib.TH1 import export1d
from lib.datacard_shapes import DataCardShapes
from modules.eftParam import EFTParam
#
import uproot
#import uproot_methods
#from ROOT import TFile, TDirectory, TH1F
#import coffea.hist as hist
import numpy as np
import pandas as pd
import json
#import matplotlib.pyplot as plt


np.random.seed(1)

'''
datacard bins defined here
use datacardshapes class to create nn bins
'''
#
#isblind = False # this is still "blinded" if true, uses NN score  up to ~.7
isblind_cat = ['partblind','partblindalt']
if len(sys.argv) > 1 :
    isblind = sys.argv[1] not in isblind_cat # # isblind false : cut out signal sensitive bins
else:
    isblind = True
#
pt_bins  = [0,200,300,450]
sdM_bins = cfg.sdm_bins     # new format
#sdM_bins = [50,75,90,105,120,140,200] # old format

#nn = 'NN'
nn = cfg.nn

#[50, 75, 85, 100, 110, 120, 135, 150, 200]
#hist2d = DataCardShapes(pt_bins,sdM_bins,n_NN_bins=10, isblind=isblind) # syntax np2d[year][ith_ptbin]
#
# initializes np.histogram2d function with defined bins
def get_sumw_sumw2(df, weights, year):
    #sumw,  _    = self.hist3d([df[nn].to_numpy(), df['Zh_M'].to_numpy(), self.ptclip(df['Zh_pt'].to_numpy())], # clip at 500 for large pt
    #                          weights=weights)
    #sumw2, _    = self.hist3d([df[nn].to_numpy(), df['Zh_M'].to_numpy(), self.ptclip(df['Zh_pt'].to_numpy())], # clip at 500 for large pt
    # weights=np.power(weights,2))
    ret_x_y = (lambda df: (df[nn].to_numpy(), df['Zh_M'].to_numpy()))
    sumw = np.array([
        hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
                             weights=weights[df['pt_bin']==i_bin])[0] 
        for i_bin in range(1,len(pt_bins))
    ])
    sumw2 = np.array([
        hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
                             weights=np.power(weights[df['pt_bin']==i_bin],2))[0]
        for i_bin in range(1,len(pt_bins))
    ])

    return sumw, sumw2
#


class MakeDataCard:
    '''
    Class to handle the overhead of creating,
    and writing the datacard as well as 
    the root file with shapes for template fit
    '''
    file_dir = './files/'
    #
    pt_bins = [0,200,300,450] # [200,300,400]
    
    dc_dir = 'Higgs-Combine-Tool'
    eft_out_file = 'EFT_Parameterization_v4.npy' # update this if settings are changed
    
    tag = '' if len(sys.argv) < 2 else sys.argv[1]+'_'
    #sig and bkg variables
    weights = ['weight','genWeight','topptWeight', #'SAT_HEMVetoWeight_drLeptonCleaned',
               'dak8md_bbvl_sf',
               #'lep_sf',
               #'lep_trigeffsf',
               'electron_sf', 'muon_sf',
               'electron_trigeffsf', 'muon_trigeffsf',
               'BTagWeight','puWeight']#,'PrefireWeight']
    weight_sys = [#'Stop0l_topptWeight_Up' ,'Stop0l_topptWeight_Down', # only for top_powheg samples
        'BTagWeight_up_jes',     'BTagWeight_down_jes',          
        'BTagWeight_up_lf',      'BTagWeight_down_lf',            
        'BTagWeight_up_hf',      'BTagWeight_down_hf',            
        'BTagWeight_up_lfstats1','BTagWeight_down_lfstats1',# non-corr
        'BTagWeight_up_lfstats2','BTagWeight_down_lfstats2',# non-corr
        'BTagWeight_up_hfstats1','BTagWeight_down_hfstats1',# non-corr
        'BTagWeight_up_hfstats2','BTagWeight_down_hfstats2',# non-corr
        'BTagWeight_up_cferr1','BTagWeight_down_cferr1',# non-corr
        'BTagWeight_up_cferr2','BTagWeight_down_cferr2',# non-corr
        #'puWeight_Up','puWeight_Down',
        'puWeightUp','puWeightDown',
        'pdfWeight_Up','pdfWeight_Down',
        'ISR_Up', 'ISR_Down',
        'FSR_Up','FSR_Down',
        'mu_r_Up','mu_r_Down',
        'mu_f_Up','mu_f_Down',
        'mu_rf_Up','mu_rf_Down',
        'topptWeight_Up' ,'topptWeight_Down',
        #'PrefireWeight_Up','PrefireWeight_Down',
        #'lep_trigeffsf_Up', 'lep_trigeffsf_Down',
        'electron_sf_up','electron_sf_down',
        'muon_sf_up','muon_sf_down',
        'electron_trigeffsf_Up', 'electron_trigeffsf_Down',
        'muon_trigeffsf_Up', 'muon_trigeffsf_Down',
        'dak8md_bbvl_sf_Up', 'dak8md_bbvl_sf_Down',

    ]
        
    bkg_v  = weights + weight_sys + [nn,'Zh_M','Zh_pt','process','sample','n_ak4jets']
    sig_v  = bkg_v + ['genZHpt']
    data_v = [nn,'Zh_M','Zh_pt','process','n_ak4jets']
    #
    if 'recoeft' not in tag:
        accepted_sig  = [f'{s}{i}' for i in range(len(pt_bins)) for s in ['ttH','ttZ']]#['ttHbb','ttZbb']] # change new and not new ttzbb here
    else:
        accepted_sig = ['ttH','ttZ']
    #accepted_bkg  = ['ttX','TTBar','tt_bb','tt_2b','VJets','other','single_t']
    #accepted_bkg  = ['ttX','TTBar','tt_bb','tt_2b','VJets','single_t']
    accepted_bkg  = ['ttX','TTBar','tt_B','VJets','single_t']
    accepted_data = ['data_obs']
    
    
    def __init__(self, 
                 sig = cfg.Sig_MC+cfg.sig_sys_samples, # new sig/bkg names
                 bkg = cfg.Bkg_MC+cfg.bkg_sys_samples, # new sig/bkg names
                 years = cfg.Years, 
                 #years = ['2017'], 
                 isblind=True,
                 sumw_sumw2=get_sumw_sumw2):
        # start getting signal and bkg 
        self.sig = sig
        self.bkg = bkg
        self.data = cfg.Data_samples
        self.years =  years
        #self.years = ['2018']
        self.isblind = isblind
        self.dc_bins = len(self.pt_bins[1:])
        self.get_sumw_sumw2 = sumw_sumw2
        #
        #self.getdata()

    def makeDC(self):
        self.getdatav2()
        if len(self.accepted_sig) > 2:
            self.process_sig(self.data_dict) # add new_ttZbb meging , need to rewrite a bit
        #print(self.data_dict.keys())
        self.initialize_hists() # switch to np.histogramdd , also, might be nice to get systematic shape plots here
        self.initialize_roofile() # format in file : $CHANNEL_$PROCESS, $CHANNEL_$PROCESS_$SYSTEMATIC
        self.initialize_datacard()
        # add systematics to histos
        self.setup_Systematics()
        self.add_Systematics()
        #
        # from lib.datacard_uncorrsys import add_uncorrSystematics
        # add_uncorrSystematics(self)
        # make eft param file is needed
        #self.make_eftparam_file()
        #self.add_wcs() # antiquated, adds eft params directly to datacard
        #
        #if self.isblind:   not necessary
        #    self.data_to_pdata()
        self.fill_roofile()
        self.close_roofile()
        self.close_dc()
        

    @t2Run
    def getdatav2(self):

        self.data_dict = {}
        #        for y in self.years:
        self.sig  = [f'{sig}_{y}'   for sig  in self.sig  for y in self.years]
        self.bkg  = [f'{bkg}_{y}'   for bkg  in self.bkg  for y in self.years]
        self.data = [f'{data}_{y}'  for data in self.data for y in self.years]
        pool = Pool(8)
        print('Starting signal worker')
        all_samples = self.sig+self.bkg+self.data
        results = pool.map(self.worker, all_samples) # to run this in parallel 
        pool.close()
        #pool.join()

        for result in results:
            if result is None: continue
            for key, value in result.items(): # should just be one thing 
                if key in self.data_dict:
                    self.data_dict[key] = pd.concat([self.data_dict[key],value], axis='rows', ignore_index=True)
                else:
                    self.data_dict[key] = value
        del results, pool
        print('Done with signal')
        #print(self.data_dict['new_tt_2b_2017'])
        #self.process_bkg(data_dict,y) # for tt_bb merging 
            
    def worker(self, process):
        y = re.search(r'201\d',process).group()
        p_vars = None # i think i have to do this to flush the memory
        p_vars = self.sig_v if process in self.sig else self.bkg_v
        p_vars = p_vars + cfg.ana_vars[f'sysvars_{y}']
        if '2018'  in process : p_vars = p_vars + ['HEM_weight']
        if 'TTBar' in process : p_vars = p_vars + ['tt_type','ttCxsecWeight','ttCxsecWeight_Up','ttCxsecWeight_Down']
        if 'ttbb' in process  : p_vars = p_vars + ['tt_type','tt2bxsecWeight','tt2bxsecWeight_Up','tt2bxsecWeight_Down']
        if 'Data'  in process : p_vars = self.data_v #[nn,'Zh_M','Zh_pt','process']
        # get important info for signal
        return self.updatedict(process.replace(f'_{y}', ''), p_vars, y)
        
    #
        
    def updatedict(self, p, v, y=''):
        sub_f_dir = 'data_files' if 'Data' in p else 'mc_files'
        if not os.path.exists(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl'): return 
        if p in cfg.all_sys_samples:
            v = [var for var in v if var not in self.weight_sys] # save memory
            df = pd.read_pickle(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl').filter(items=v)
        else:
            df = pd.read_pickle(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl').filter(items=v)
            if 'Data' not in p: # clip mu_rf, isr/fsr, pdf at 3sigma percentile, 99.7% (.15%,99.85%)
                func = np.nanpercentile # because this is a long name
                for ud in ['Up','Down']:
                    for v_str in [ f'{s}_{ud}' for s in ['mu_r','mu_f','mu_rf','ISR','FSR','pdfWeight']]:
                        df.loc[:,v_str] = df[v_str].clip(func(df[v_str].values,.15), func(df[v_str].values,99.85))
        df.loc[:,'Zh_pt'] = df['Zh_pt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]+1)
        df['pt_bin'] = pd.cut(df['Zh_pt'], bins=self.pt_bins+[500],
                              labels=[i_bin for i_bin in range(len(self.pt_bins))])
        #df = df[df[nn] >= 0.0]
        df = df[(df[nn] >= 0.0) & (df['n_ak4jets'] >=5)] # trying >= 5 cut
        group = df.groupby(by='process')
        if 'Data' not in p:
            for n,g in group:
                if n not in ['ttH','ttZ'] + self.accepted_bkg: 
                    df.drop(g.index, inplace=True)
        group = df.groupby(by='process')
        del df
        # extract sys type (if any)
        sys = '' 
        data_dict = {}
        if p in cfg.all_sys_samples: 
            if __name__ != '__main__': return data_dict
            sys =  '_'+p.split('_')[-1] # format [process_name]_[systype]
            if 'sys' in p: # hdamp and UE here
                data_dict.update({f"{n}_{y}_hdamp{'_ttbb' if 'ttbb' in p else ''}Up": g[g['sample'].str.contains('hdampUp')] for n,g in group})
                data_dict.update({f"{n}_{y}_hdamp{'_ttbb' if 'ttbb' in p else ''}Down": g[g['sample'].str.contains('hdampDown')] for n,g in group})
                if 'ttbb' not in p:
                    data_dict.update({f'{n}_{y}_UEUp': g[g['sample'].str.contains('UEUp')] for n,g in group})
                    data_dict.update({f'{n}_{y}_UEDown': g[g['sample'].str.contains('UEDown')] for n,g in group})
                return data_dict
            #if 'h' in sys and 'ttbb' in p: sys = sys.replace('hdamp','hdamp_ttbb') 
            if 'jmr' in sys or 'jms' in sys or 'jer' in sys :   
                sys = sys.replace('Up',f'_{y}Up').replace('Down',f'_{y}Down')
            
        data_dict = {f"{n.replace('Data','data_obs')}_{y}{sys}": g for n,g in group} # iterate over name and content
        return data_dict

    def process_sig(self, data_dict, y=''):
        sig_groups = ['ttZ', 'ttH']#, 'new_ttZbb'] 
        sig_names = []
        for g in sig_groups:
            sig_names = sig_names + re.findall(rf'{g}_201\d\w*', ' '.join(data_dict.keys()))
        # should just use findall to get all relevent samples to break up by genZHpt
        for sig_name in sig_names:
            #sig_name = f'{s}_{y}'
            s_pre = re.search(r'(ttH|ttZ)',sig_name).group()
            #print(sig_name)
            if sig_name not in data_dict: continue
            df = data_dict[sig_name]
            for i,_ in enumerate(self.pt_bins[:-1]):
                # insert 1 after sig_goup name # format should be [siggroup]_[sys]_[year]
                new_sig_name = f"{s_pre}{i}{sig_name.replace(s_pre, '')}"
                #print(new_sig_name)
                data_dict[new_sig_name] = df[(df['genZHpt'] >= self.pt_bins[i]) & 
                                             (df['genZHpt'] < self.pt_bins[i+1])]
                #
            #data_dict[f'{s}{len(pt_bins) -1 }_{y}'] = df[df['genZHpt'] >= pt_bins[-1]]
            data_dict[f"{s_pre}{len(self.pt_bins) -1}{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= self.pt_bins[-1]]
            data_dict.pop(sig_name)
    
    def process_bkg(self, data_dict, y): 
        pass # this will eventually be the merging procedure for new_tt_bb
    
    @t2Run
    def initialize_hists(self):
        # create 3d histogram
        # datacard is binned in pt
        # and tamplate root file is a "flattened"
        # 2d histo of Zh_M and NN score
        # initialize and 'fill' here
        #
        '''
        Well, this might be a bit too simplistic, im thinking create a 2d histogram dictionary, per year and per reco pt channel
        histo2d_dict = {'2016':{'pt1': 2d_hist_object}}
        '''
        #self.hist3d = functools.partial(
        #    np.histogramdd, 
        #    bins=[self.nn_bins,[50,80,105,145,200],pt_bins[1:]+[500]]) ## skip 0-200 bin since it should be empty add 500 to account for >= pt_bins[-1]
        
        #self.ptclip = functools.partial(np.clip, a_min=pt_bins[1]+1, a_max=499.)
        self.histos = {}
        edges = []
        for s,v in self.data_dict.items():
            #y = s.split('_')[-1] # format blah_blah_"year" 
            #print(s)
            y = re.search(r'201\d',s).group() # better than previous method
            w = getZhbbWeight(v,y) if 'data' not in s else np.ones_like((v[nn] if nn in v else v.iloc[:,0]).to_numpy()) # which should be 1 to handle data
            sumw, sumw2 = self.get_sumw_sumw2(v, w, y)
            self.histos[s] = {'sumw':sumw, # dim 1,2,3 is zh_pt, nn, zh_m
                              'sumw2':sumw2}
        

    def setup_Systematics(self):
        process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        Systematic.set_dc_processes(self.dc_dict, process_line)
        ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos)#, self.hist3d, self.ptclip)

    @t2Run
    def add_Systematics(self):
        # SIG = ttZbb[0,1,2,3] ttHbb[0,1,2,3]
        # BKG = ttX, TTBar, old_tt_bb, VJets, other
        #process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        all_mc = self.accepted_sig + self.accepted_bkg
        all_but_ttbb = self.accepted_sig + ['TTBar','ttX','VJets','single_t']
        #tth_sig  = re.findall(r'ttH\d', ' '.join(self.accepted_sig))
        #ttz_sig  = re.findall(r'\w*ttZ\d', ' '.join(self.accepted_sig))
        tth_sig = [s for s in self.accepted_sig if 'ttH' in s]
        ttz_sig = [s for s in self.accepted_sig if 'ttZ' in s]
        ttbar_mc = ['TTBar','tt_B']
        jec_mc   = ttbar_mc + ['single_t'] + tth_sig + ttz_sig
        #
        #Systematic.set_dc_processes(self.dc_dict, process_line)
        self.write2dc(f'# Rate uncertainties\n')
        # new lumi
        Systematic('lumi_13TeV_2016',       'lnN',  all_mc, 1.01)
        Systematic('lumi_13TeV_2017',       'lnN',  all_mc, 1.02)
        Systematic('lumi_13TeV_2018',       'lnN',  all_mc, 1.015)
        Systematic('lumi_13TeV_correlated', 'lnN',  all_mc, {'2016':1.006, '2017': 1.009, '2018': 1.02}) # 0.6, 0.9, 2.0
        Systematic('lumi_13TeV_1718',       'lnN',  all_mc, {'2017':1.006, '2018': 1.002}) # 0.6, 0.2
        # old 
        #Systematic('lumi_2016', 'lnN',  all_mc, 1.025)
        #Systematic('lumi_2017', 'lnN',  all_mc, 1.023)
        #Systematic('lumi_2018', 'lnN',  all_mc, 1.025)
        
        # probably do xsec theo rates here
        # first signal
        Systematic('tth_ggpdf', 'lnN', tth_sig, 1.036)       # revisit these
        Systematic('ttz_ggpdf', 'lnN', ttz_sig, 1.035)       # revisit these
        Systematic('tth_qsc' ,  'lnN', tth_sig, 1.058,0.908) # revisit these
        Systematic('ttz_qsc'  , 'lnN', ttz_sig, 1.081,0.907) # revisit these
        # now background                                     # revisit these
        Systematic('ggpdf', 'lnN', ttbar_mc, 1.042)          # revisit these
        #Systematic('qqpdf', 'lnN', ['ttX','VJets','other'],
        #           [p_norms.rate_unc['pdf']['ttX'],
        #            p_norms.rate_unc['pdf']['VJets'],
        #            p_norms.rate_unc['pdf']['other']])# combine into one 
        Systematic('qqpdf', 'lnN', ['ttX','VJets'],
                   [p_norms.rate_unc['pdf']['ttX'],
                    p_norms.rate_unc['pdf']['VJets']])
        Systematic('qgpdf', 'lnN', ['single_t'], [p_norms.rate_unc['pdf']['single_t']])
        #Systematic('qqpdf', 'lnN', ['VJets'], 1.038)# combine into one
        #Systematic('qqpdf', 'lnN', ['other'], 1.050)# combine into one
        #Systematic('tt_qsc'   , 'lnN', ttbar_mc+['ttX'], [[1.024,0.965] for _ in ttbar_mc]+[1.300])
        Systematic('tt_qsc'   , 'lnN', ttbar_mc, [[1.024,0.965] for _ in ttbar_mc])
        Systematic('ttx_qsc'  ,  'lnN', ['ttX'], p_norms.rate_unc['QCD_scale']['ttX']) 
        Systematic('singlet_qsc'  ,  'lnN', ['single_t'], p_norms.rate_unc['QCD_scale']['single_t']) 
        Systematic('v_qsc'    , 'lnN', ['VJets'],p_norms.rate_unc['QCD_scale']['VJets'])#1.008, 0.996) # .821/1.24
        #Systematic('other_qsc', 'lnN', ['other'],p_norms.rate_unc['QCD_scale']['other']) #1.05, 0.95)   # .898/1.12
        #Systematic('tt2bxsec', 'lnN', ['old_tt_bb'], 1.5)
        # Shape Systatics
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Shape uncertainties \n')
        #ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos)#, self.hist3d, self.ptclip)
        # when defining shape, must incluse whether it is a mcsta, scale, or up/down syst
        for y in self.years:
            #self.histos = ShapeSystematic(f'btglf_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_lf', 'BTagWeight_down_lf').get_shape()
            #self.histos = ShapeSystematic(f'btghf_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_hf', 'BTagWeight_down_hf').get_shape()
            self.histos = ShapeSystematic(f'btglfstats1_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_lfstats1','BTagWeight_down_lfstats1').get_shape()
            self.histos = ShapeSystematic(f'btglfstats2_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_lfstats2','BTagWeight_down_lfstats2').get_shape()
            self.histos = ShapeSystematic(f'btghfstats1_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_hfstats1','BTagWeight_down_hfstats1').get_shape()
            self.histos = ShapeSystematic(f'btghfstats2_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_hfstats2','BTagWeight_down_hfstats2').get_shape()
            #self.histos = ShapeSystematic(f'lepsf_{y}', 'shape', 'up/down', all_mc, 1, 'lep_sf_up','lep_sf_down').get_shape()
            self.histos = ShapeSystematic(f'eleclepsf_{y}', 'shape', 'up/down', all_mc, 1, 'electron_sf_up','electron_sf_down').get_shape()
            self.histos = ShapeSystematic(f'mulepsf_{y}', 'shape', 'up/down', all_mc, 1, 'muon_sf_up','muon_sf_down').get_shape()
            #self.histos = ShapeSystematic(f'trigeffsf_{y}', 'shape', 'up/down', all_mc, 1, 'lep_trigeffsf_Up', 'lep_trigeffsf_Down').get_shape()
            self.histos = ShapeSystematic(f'electrigeffsf_{y}', 'shape', 'up/down', all_mc, 1, 'electron_trigeffsf_Up', 'electron_trigeffsf_Down').get_shape()
            self.histos = ShapeSystematic(f'mutrigeffsf_{y}', 'shape', 'up/down', all_mc, 1, 'muon_trigeffsf_Up', 'muon_trigeffsf_Down').get_shape()
            self.histos = ShapeSystematic(f'pu_{y}',      'shape', 'up/down', all_mc, 1, 'puWeightUp','puWeightDown', extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'bbvlsf_{y}',      'shape', 'up/down', all_mc, 1, 'dak8md_bbvl_sf_Up', 'dak8md_bbvl_sf_Down').get_shape()            
            # new jes and jer
            self.histos = ShapeSystematic(f'jesRelativeSample{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'jesHF{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'jesAbsolute{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'jesEC2{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'jesBBEC1{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'jms_{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'jmr_{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            #self.histos = ShapeSystematic(f'jmssig_{y}', 'shape', 'qconly', tth_sig+ttz_sig, 1, extraQC=True).get_shape()
            #self.histos = ShapeSystematic(f'jmrsig_{y}', 'shape', 'qconly', tth_sig+ttz_sig, 1, extraQC=True).get_shape()
            #self.histos = ShapeSystematic(f'jmsttbb_{y}', 'shape', 'qconly', ['tt_B'], 1, extraQC=True).get_shape()
            #self.histos = ShapeSystematic(f'jmrttbb_{y}', 'shape', 'qconly', ['tt_B'], 1, extraQC=True).get_shape()
            #self.histos = ShapeSystematic(f'jmstt_{y}', 'shape', 'qconly', ['single_t','TTBar'], 1, extraQC=True).get_shape()
            #self.histos = ShapeSystematic(f'jmrtt_{y}', 'shape', 'qconly', ['single_t','TTBar'], 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'ak4jer_{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            self.histos = ShapeSystematic(f'ak8jer_{y}', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        # new jes and jer
        self.histos = ShapeSystematic('jesRelativeBal', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('jesFlavorQCD', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('jesHF', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('jesAbsolute', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('jesEC2', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('jesBBEC1', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('jesHEMIssue', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
        #
        self.histos = ShapeSystematic(f'btglf', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_lf', 'BTagWeight_down_lf').get_shape()
        self.histos = ShapeSystematic(f'btghf', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_hf', 'BTagWeight_down_hf').get_shape()
        self.histos = ShapeSystematic(f'btgjes', 'shape', 'up/down', all_mc, 1,'BTagWeight_up_jes','BTagWeight_down_jes').get_shape()
        self.histos = ShapeSystematic(f'btgcferr1', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_cferr1','BTagWeight_down_cferr1').get_shape()
        self.histos = ShapeSystematic(f'btgcferr2', 'shape', 'up/down', all_mc, 1, 'BTagWeight_up_cferr2','BTagWeight_down_cferr2').get_shape()
        #self.histos = ShapeSystematic(f'lepsf', 'shape', 'up/down', all_mc, 1, 'lep_sf_up','lep_sf_down').get_shape()
        self.histos = ShapeSystematic(f'pref_2016', 'shape', 'up/down', all_mc, 1, 'PrefireWeight_Up' ,'PrefireWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pref_2017', 'shape', 'up/down', all_mc, 1, 'PrefireWeight_Up' ,'PrefireWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'toppt', 'shape', 'up/down', ttbar_mc, 1, 'topptWeight_Up' ,'topptWeight_Down').get_shape() # using hacky unc.
        self.histos = ShapeSystematic(f'isr_tt', 'shape', 'ps', ['TTBar'], 1, 'ISR_Up','ISR_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'fsr_tt', 'shape', 'ps', ['TTBar'], 1, 'FSR_Up','FSR_Down', extraQC=True).get_shape()
        #self.histos = ShapeSystematic(f'mu_r', 'shape', 'normup/down', all_but_ttbb, 1, 'mu_r_Up','mu_r_Down').get_shape()
        #self.histos = ShapeSystematic(f'mu_f', 'shape', 'normup/down', all_but_ttbb, 1, 'mu_f_Up','mu_f_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_r_tt', 'shape', 'normup/down', ['TTBar'], 1, 'mu_r_Up','mu_r_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_f_tt', 'shape', 'normup/down', ['TTBar'], 1, 'mu_f_Up','mu_f_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_r_tth', 'shape', 'normup/down', tth_sig, 1, 'mu_r_Up','mu_r_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'mu_f_tth', 'shape', 'normup/down', tth_sig, 1, 'mu_f_Up','mu_f_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'mu_r_ttz', 'shape', 'normup/down', ttz_sig, 1, 'mu_r_Up','mu_r_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'mu_f_ttz', 'shape', 'normup/down', ttz_sig, 1, 'mu_f_Up','mu_f_Down', extraQC=True).get_shape()
        #
        self.histos = ShapeSystematic(f'isr_tth', 'shape', 'normup/down', tth_sig, 1, 'ISR_Up','ISR_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'isr_ttz', 'shape', 'normup/down', ttz_sig, 1, 'ISR_Up','ISR_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'fsr_tth', 'shape', 'normup/down', tth_sig, 1, 'FSR_Up','FSR_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic(f'fsr_ttz', 'shape', 'normup/down', ttz_sig, 1, 'FSR_Up','FSR_Down', extraQC=True).get_shape()
        #self.histos = ShapeSystematic(f'isr_ttbb', 'shape', 'ps', ['tt_bb', 'tt_2b'], 1, 'ISR_Up','ISR_Down').get_shape()
        #self.histos = ShapeSystematic(f'fsr_ttbb', 'shape', 'ps', ['tt_bb', 'tt_2b'], 1, 'FSR_Up','FSR_Down', extraQC=True).get_shape()
        #self.histos = ShapeSystematic(f'mu_r_ttbb', 'shape', 'normup/down', ['tt_bb', 'tt_2b'], 1, 'mu_r_Up','mu_r_Down').get_shape()
        #self.histos = ShapeSystematic(f'mu_f_ttbb', 'shape', 'normup/down', ['tt_bb', 'tt_2b'], 1, 'mu_f_Up','mu_f_Down').get_shape()
        self.histos = ShapeSystematic('isr_ttbb', 'shape', 'ps', ['tt_B'], 1, 'ISR_Up','ISR_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic('fsr_ttbb', 'shape', 'ps', ['tt_B'], 1, 'FSR_Up','FSR_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic('mu_r_ttbb', 'shape', 'normup/down', ['tt_B'], 1, 'mu_r_Up','mu_r_Down', extraQC=True).get_shape()
        self.histos = ShapeSystematic('mu_f_ttbb', 'shape', 'normup/down', ['tt_B'], 1, 'mu_f_Up','mu_f_Down', extraQC=True).get_shape()
        # redundant #self.histos = ShapeSystematic(f'mu_rf', 'shape', 'normup/down', all_mc, 1, 'mu_rf_Up','mu_rf_Down').get_shape()
        #
        #self.histos = ShapeSystematic(f'pdf_ttz', 'shape', 'up/down', ttz_sig, 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        #self.histos = ShapeSystematic(f'pdf', 'shape', 'normup/down', all_but_ttbb, 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pdf', 'shape', 'normup/down', tth_sig+ttz_sig+['TTBar'], 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        #self.histos = ShapeSystematic(f'pdf_ttbb', 'shape', 'normup/down', ['tt_bb','tt_2b'],                    1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pdf_ttbb', 'shape', 'normup/down', ['tt_B'],                    1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        
        # These shapes are already computed, just need to add to datacard
        self.histos = ShapeSystematic('UE',     'shape', 'qconly', ['TTBar'], 1, extraQC=True).get_shape()
        #self.histos = ShapeSystematic('erdOn', 'shape',  ['TTBar'], 1) # not working with just one shape at the moment
        self.histos = ShapeSystematic('hdamp', 'shape', 'qconly', ['TTBar'], 1, extraQC=True).get_shape()
        #self.histos = ShapeSystematic('hdamp_ttbb', 'shape', 'qconly', ['tt_bb', 'tt_2b'], 1, extraQC=True).get_shape()
        self.histos = ShapeSystematic('hdamp_ttbb', 'shape', 'qconly', ['tt_B'], 1, extraQC=True).get_shape()
        #
        self.write2dc(100*'-'+'\n')
        self.histos = ShapeSystematic('tt2bxsec', 'shape', 'up/down', ['tt_B'],  1, 'tt2bxsecWeight_Up', 'tt2bxsecWeight_Down').get_shape()
        self.histos = ShapeSystematic('ttCxsec',  'shape', 'up/down', ['TTBar'], 1, 'ttCxsecWeight_Up',  'ttCxsecWeight_Down').get_shape()
        self.write2dc('# Float tt_B normalization\n') 
        self.write2dc('CMS_ttbbnorm rateParam * tt_B 1 [0.0,5.0]\n')
        #Systematic('CMS_ttbbnorm', 'lnN', ['tt_bb','tt_2b'], 10)
        self.write2dc(100*'-'+'\n')
        self.write2dc('# MC Stats uncertainties\n') 
        #self.histos = ShapeSystematic(f'mcstat','shape','mcstat', all_mc, 1).get_shape()
        self.write2dc('* autoMCStats 10 0  1\n') 
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Group definitions \n') 
        self.write2dc('sig_rtheo group = tth_qsc ttz_qsc tth_ggpdf ttz_ggpdf\n')
        self.write2dc('theo group = CMS_ttbbnorm tt2bxsec ttCxsec hdamp_ttbb hdamp UE toppt pdf_ttbb pdf\n')
        self.write2dc('theo group += mu_f_ttbb mu_r_ttbb mu_f_tt mu_r_tt mu_f_tth mu_r_tth mu_f_ttz mu_r_ttz\n')
        self.write2dc('theo group += isr_ttbb fsr_ttbb isr_tt fsr_tt isr_tth fsr_tth isr_ttz fsr_ttz\n')
        self.write2dc('theo group += tth_ggpdf ttz_ggpdf tth_qsc ttz_qsc ggpdf qqpdf qgpdf tt_qsc ttx_qsc singlet_qsc v_qsc\n')
        #

    def make_eftparam_file(self):
        out_path = sys.path[1]+'Higgs-Combine-Tool/'
        out_file = self.eft_out_file #'EFT_Parameterization_v2.npy' # manually change this 
        if os.path.exists(out_path+out_file): return
        eft = EFTParam()
        #eft.save_to_dict(year='2016', force_year='2016')
        #eft.save_to_dict(year='2017', force_year='2017')
        #eft.save_to_dict(year='2018', force_year='2018')
        eft.save_to_dict(year='2016', force_year='2018')
        eft.save_to_dict(year='2017', force_year='2018')
        eft.save_to_dict(year='2018', force_year='2018')
        import pickle
        pickle.dump(eft.out_dict, open(out_path+out_file,'wb'), protocol=2)

    def add_wcs(self):
        eft = EFTParam()
        for y in self.years:
            dc_lines = eft.get_EFT_lines(year='2018') # adds dc lines as an attribute to EFTParam class # add 2018 calc to all dc
            self.writelines2dc(dc_lines, y)

    def write2dc(self, str2write):
        for y in self.years:
            self.dc_dict[y].write(str2write)
    
    def writelines2dc(self,lines,y): # by year
        if lines: 
            self.dc_dict[y].writelines(lines)
            
    def data_to_pdata(self):
        for y in self.years:
            sumw  = np.nansum(np.array([self.histos[f'{p}_{y}']['sumw']  for p in self.accepted_sig + self.accepted_bkg if f'{p}_{y}' in self.histos] ), axis=0)
            #sumw2 = np.power(sumw,2)
            self.histos[f'data_obs_{y}'] = {'sumw' : sumw,
                                            'sumw2': sumw}

    @t2Run
    def fill_roofile(self):
        for p,v in self.histos.items():
            y = re.findall(r'201\d',p)[0] # first instance of this should be the process year
            process = p.split(f'_{y}')[0]
            if process not in self.accepted_sig + self.accepted_bkg + self.accepted_data: continue
            if p.replace(f'_{y}', '') not in self.accepted_sig + self.accepted_bkg + self.accepted_data:
                # this should mean its a shape systematic or a process not in accepted
                #process = p.split(f'_{y}')[0]
                sys     = p.replace(f'{process}_{y}', '') # should have format _sys
            else:
                #process = p.replace(f'_{y}' ,'')
                sys     = ''
            #
            for pt_bin in range(v['sumw'].shape[0]):
                if pt_bin == 0: # 0,1,2 (-1)
                    if not self.isblind:
                        to_flat = (lambda a: self.take_first_last_mbins(self.merge_last_mbin(pt_bin,a)))
                    else:
                        to_flat = (lambda a: self.merge_last_mbin(pt_bin,a).flatten())
                else: # not the first pt bin
                    if not self.isblind:
                        to_flat = (lambda a : self.take_first_last_mbins(a[pt_bin,:,:]))
                    else:
                        to_flat = (lambda a : a[pt_bin,:,:].flatten())
                        #to_flat = (lambda a : self.merge_lastn_nnbin(pt_bin, a))
                temp_dict = {'sumw' : to_flat(v['sumw'])}#* (1 if y != '2017' else 3.3032)}#2.2967)} # to just scale to full run2
                #temp_dict = {'sumw' : to_flat(v['sumw'])* (1 if y != '2018' else cfg.Lumi['run2']/cfg.Lumi['2018'])} # to just scale to full run2
                hist_name = f'Zhpt{pt_bin+1}_{process}{sys}'
                if 'sumw2' in v:
                    temp_dict['sumw2'] = to_flat(v['sumw2']) 
                    #temp_dict['sumw2'] = to_flat(v['sumw2']) if 'tt_B' not in p else to_flat(v['sumw2'])/4 
                self.roo_dict[y][hist_name] = export1d(temp_dict, hist_name)

    @staticmethod
    def merge_last_mbin(pt_bin,a):
        a[pt_bin,:,-2] = a[pt_bin,:,-2] + a[pt_bin,:,-1]
        #return a[pt_bin,:,:-1].flatten()
        return a[pt_bin,:,:-1]
    #@staticmethod
    #def take_first_last_mbins(pt_bin,a): # but not for the first 3 NN quantile
    #    sb = np.stack([a[pt_bin,3:,0],a[pt_bin,3:,-1]], axis=-1)
    #    return np.append(a[pt_bin,0:3,:], sb).flatten()
    @staticmethod
    def take_first_last_mbins(a): # but not for the first 3 NN quantile
        sb = np.stack([a[3:,0],a[3:,-1]], axis=-1)
        return np.append(a[0:3,:], sb).flatten()

    @staticmethod
    def merge_lastn_nnbin(pt_bin,a):
        n_ = 2 # merge last 2
        a[pt_bin,-1*n_,:] = sum(a[pt_bin,-1*n_:,:])
        return a[pt_bin,:-1*(n_-1),:].flatten()
    @staticmethod
    def merge_to_4mbin(pt_bin,a):
        print(a[pt_bin].flatten())
        a[pt_bin,7:,1] = np.sum(a[pt_bin,7:,1:3],axis=1)   # concat last set of bins
        a[pt_bin,7:,4] = np.sum(a[pt_bin,7:,4:],axis=1)   # concat last set of bins
        a = np.append(a[pt_bin,:7,:],np.hstack([a[pt_bin,7:,0:2],a[pt_bin,7:,4:]]))
        print(a)
        return a
    #@staticmethod
    #def merge_last_nnbin(pt_bin,a):
    #    print(a[pt_bin].flatten())
    #    a[pt_bin,-2,:] = sum(a[pt_bin,-2:,:])   # concat last set of bins
    #    a[pt_bin,-4,:] = sum(a[pt_bin,-4:-2,:]) # concat 2nd to last set of bins
    #    a= np.append(a[pt_bin,:-3,:], a[pt_bin,-2,:] ) # this will be flattened automatically
    #    print(a)
    #    return a

    def initialize_datacard(self):
        # creat 3 data cards , 1 per year
        dc_dict = {y: open(f'{self.dc_dir}/datacard_{self.tag}{y}.txt', 'w') for y in self.years}
        for y,txt in dc_dict.items():
            txt.writelines([
                f'Datacard for {y}\n',
                'imax * number of bins\n',
                'jmax * number of processes minus 1\n',
                'kmax * number of nuisance paramerters\n',
                100*'-'+'\n',
                f'shapes data_obs * datacard_{self.tag}{y}.root $CHANNEL_data_obs\n',
                f'shapes * * datacard_{self.tag}{y}.root $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC\n',
                100*'-'+'\n',
                f"{'bin':20}{' '.join(['Zhpt'+str(i+1) for i in range(self.dc_bins)])}\n",
                f"{'observation':20}{self.dc_bins*'-1  '}\n",
                100*'-'+'\n',
                f"{'bin':20}{' '.join(['Zhpt'+str(i+1) for i in range(self.dc_bins) for _ in range(len(self.accepted_sig + self.accepted_bkg))])}\n",
                f"{'process':20}{' '.join(s for _ in range(self.dc_bins) for s in self.accepted_sig+self.accepted_bkg)}\n",
                f"{'process':20}{' '.join(str(i) for _ in range(self.dc_bins) for i in range(-len(self.accepted_sig)+1, len(self.accepted_bkg)+1))}\n",
                f"{'rate':20}{' '.join(str(-1) for _ in range(self.dc_bins) for s in self.accepted_sig + self.accepted_bkg)}\n",
                100*'-'+'\n'])
        #
        self.dc_dict = dc_dict

    def initialize_roofile(self):
        # create root file with uproot methods
        # for this use case
        roo_dict = {}
        for y in self.years:
            roo_name = f'{self.dc_dir}/datacard_{self.tag}{y}.root'
            if os.path.exists(roo_name):
                os.system(f"rm {roo_name}")
            roo_dict[y] = uproot.create(roo_name) 
        self.roo_dict = roo_dict

    def close_roofile(self):
        for roo in self.roo_dict:
            self.roo_dict[roo].close()

    def close_dc(self):
        for dc in self.dc_dict:
            self.dc_dict[dc].close()

class Systematic: # Class to handle Datacard systematics 
    ''' 
    Syntax: Systematic(Systematic Name,Systematic Type, Channel, Affected Processes, Value, Additional Information)
    '''
    #dc_root_dir = 'Higgs-Combine-Tool/'
    dc_root_dir = ''
    datacard     = None
    allowed_processes = None
    

    def __init__(self, name, stype, process_ids, value, optvalue=None, info=None):
        self.name     = name
        self.stype    = stype
        self.ids      = process_ids
        #self.channel  = channel
        self.years     = re.findall(r'201\d', name) if name != 'jesHEMIssue' else ['2018']
        if len(self.years) == 0: self.years = cfg.Years
        # handle different value casts
        if type(value) is list:
            self.value = {i:v for i,v in zip(self.ids,value)}
        elif type(value) is float or type(value) is int:
            self.value = value
        elif type(value) is dict:
            self.value = None
            self.years = list(value.keys())
        else:
            self.value = None
        #self.value    = value if type(value) is not list else {i:v for i,v in zip(self.ids,value)}
        self.optvalue = optvalue
        self.info     = '' if info is None else info
        #
        if self.datacard is not None: 
            for year in self.years:
                if type(value) is dict:
                    self.value = value[year]
                self.datacard[year].write(self.get_DC_line()) # write to datacard file upon instance creation

    @property
    def line(self):
        return '{0:14} {1:6}'.format(self.name,self.stype)

    @classmethod
    def set_dc_processes(cls,datacard, processes):
        cls.datacard = datacard
        cls.allowed_processes = processes
        cls.p_norms = json.load(open(f'{cfg.dataDir}/process_norms/process_norms_ttbbw_run2.json','r'))

    def get_DC_line(self):
        _line = self.line
        value = self.value
        optvalue = self.optvalue
        for p in self.allowed_processes:
            #_process = p.replace('\t', '').replace(' ','' )# reformat process to exclude \t 
            if p in self.ids: 
                if type(self.value) is dict:
                    value = self.value[p]
                    if type(value) is list:
                        value, optvalue = value[0], value[1]
                    else:
                        value, optvalue = value, None
                if optvalue is None:
                    _line +='{0:12}'.format(str(value))
                else:
                    entry = '{1}/{0}'.format(str(value),str(optvalue))
                    _line += f'{entry:12}'
            else :
                _line +='{0:12}'.format('-')
        _line += '\t'+self.info+'\n'
        return _line


class ShapeSystematic(Systematic): # Class to handle Datacard shape systematics
    
    '''
    Create additional histograms with systematic variations and adds them to histo dict
    Supports MCStats
    Inheirits from the Systematics Class
    Need to set cut_op and bins_dict to make effective use of this class
    '''

    bins_dict = None
    cut_op   = None
    df       = None

    def __init__(self, name, stype, subtype, process_ids, value, up=None, down=None, info=None, extraQC=False):
        super().__init__(name, stype, process_ids, value, info)
        self.up     = up
        self.down   = down
        self.subtype = subtype
        #self.extraQC = False # to make smoothing before and after plots
        self.extraQC = extraQC

    def get_shape(self):
        fun_dict = {'mcstat' :self.makeMCStatHist,
                    #'scale'  :self.makeScaleHist,
                    'up/down'    :self.makeUpDownHist,
                    'normup/down':self.makeUpDownHist,
                    'ps'     :self.makeUpDownHist,#self.handlePSUpDown,
                    'erdOn'  :self.handleErdOn} # unique case where PS weights for 2017 tth and PS weight for 2016 and 2017 need to be handled 
        if self.subtype in fun_dict: fun_dict[self.subtype]()
        #
        #if False: # testing qc vs no qc 
        if self.extraQC and __name__ == '__main__': self.handle_extraQC()
        self.handleQC()
        #
        return self.histos

    @classmethod
    def set_df_histos_histfuncs(cls, data_dict, hist_dict, 
                                sumw_sumw2=get_sumw_sumw2):# hist_func, ptclip):
        cls.data_dict = data_dict
        cls.histos = hist_dict # to basically add histos dict to scope
        cls.get_sumw_sumw2 = staticmethod(sumw_sumw2)


    def handleErdOn(self):
        for process in self.ids:
            for y in self.years:
                df = self.data_dict[f'{process}_{y}_erdOn']
                weight =  getZhbbWeight(df,y)
                sumw, sumw2 = get_sumw_sumw2(df,w_var,y)
                self.histos[f'{process}_{y}_{self.name}'] = {'sumw':sumw, 
                                                             'sumw2': sumw2}
                

    def handlePSUpDown(self):
        self.makeUpDownHist() # add histos to histos dict
        # special handling for ttHbb[0,1,2,3]_2017 and ttZbb[0,1,2,3]_[2016,2017]
        for ud_str, ud in zip(['Up','Down'],[self.up,self.down]):
            for i in range(4):
                sumw_var = ((self.histos[f'ttHbb{i}_2016_{self.name}{ud_str}']['sumw']/self.histos[f'ttHbb{i}_2016']['sumw']) +  # we want the average deviation
                            (self.histos[f'ttHbb{i}_2018_{self.name}{ud_str}']['sumw']/self.histos[f'ttHbb{i}_2018']['sumw']))/2
                sumw2_var = ((self.histos[f'ttHbb{i}_2016_{self.name}{ud_str}']['sumw2']/self.histos[f'ttHbb{i}_2016']['sumw2']) + 
                            (self.histos[f'ttHbb{i}_2018_{self.name}{ud_str}']['sumw2']/self.histos[f'ttHbb{i}_2018']['sumw2']))/2
                self.histos[f'ttHbb{i}_2017_{self.name}{ud_str}'] = {'sumw': np.nan_to_num(sumw_var) * self.histos[f'ttHbb{i}_2017']['sumw'],
                                                                     'sumw2': np.nan_to_num(sumw2_var) * self.histos[f'ttHbb{i}_2017']['sumw2']}

                #sumw_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw']/np.where(self.histos[f'ttZbb{i}_2018']['sumw'] <= 0, 
                #0.00001,self.histos[f'ttZbb{i}_2018']['sumw'])
                sumw_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw']/self.histos[f'ttZbb{i}_2018']['sumw']

                sumw2_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw2']/self.histos[f'ttZbb{i}_2018']['sumw2']
                #sumw2_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw2']/np.where(self.histos[f'ttZbb{i}_2018']['sumw2']<=0, 
                #0.00001, self.histos[f'ttZbb{i}_2018']['sumw2'])
                handles_div0 = (lambda x: np.nan_to_num(np.where( ((abs(x)>3.0) | (abs(self.histos[f'ttZbb{i}_2018']['sumw']) <0.05)) ,1, x)))
                self.histos[f'ttZbb{i}_2016_{self.name}{ud_str}'] = {'sumw'  :  handles_div0(sumw_var) * self.histos[f'ttZbb{i}_2016']['sumw'],
                                                                     'sumw2' :  handles_div0(sumw2_var) * self.histos[f'ttZbb{i}_2016']['sumw2']}
                self.histos[f'ttZbb{i}_2017_{self.name}{ud_str}'] = {'sumw'  :  handles_div0(sumw_var) * self.histos[f'ttZbb{i}_2017']['sumw'],
                                                                     'sumw2' :  handles_div0(sumw2_var) * self.histos[f'ttZbb{i}_2017']['sumw2']}
                
    def makeUpDownHist(self):
        for process in self.ids:
            for y in self.years:
                if f'{process}_{y}' not in self.data_dict: continue
                df = self.data_dict[f'{process}_{y}']
                #w_nom  = self.up.replace('_'+self.up.split('_')[-1],'') # format should be weightName_up/down
                w_nom  = self.up.replace(re.search(r'(Up|up)\w*',self.up).group(),'').rstrip('_') # format should be weightName_up/down_foobaropt
                #print(w_nom)
                #print(self.name, process)
                if w_nom == 'pdfWeight' or 'mu_' in w_nom or self.subtype == 'ps' or 'ISR' in w_nom or 'FSR' in w_nom: 
                    df[w_nom] = 1.0
                weight =  getZhbbWeight(df,y)
                for ud_str, ud in zip(['Up','Down'],[self.up,self.down]):
                    nominal_weight = weight
                    #if 'new_tt_' in process and self.subtype is 'ps':
                        #nominal_weight = (weight/df['weight']) * df[f'{ud}_weight'] # to insert correct event weight
                    #    if 'mu_' in w_nom: df[w_nom] = df[f'{ud}_weight'] # to fix issue with tt_bb
                    w_var = nominal_weight*df[ud]/df[w_nom]
                    if 'norm' in self.subtype or ('tt_' in process and self.subtype == 'ps'): 
                        #p_norms_key = re.sub(r'tt_','tt', process) if 'tt_bb' in process or 'tt_2b' in process else re.sub(r'\d$','', process)
                        p_norms_key = re.sub(r'tt_B','ttbb', process) if 'tt_B' in process else process#re.sub(r'\d$','', process)
                        w_var = w_var * self.p_norms[y][p_norms_key][ud]#sum(weight)/sum(w_var) # to keep the nominal normalization
                    #
                    #if 'new_tt_' in process and ('mu_' in w_nom or self.subtype is 'ps'):
                    #    print(process, y, ud, ud_str)
                    #    print(sum(weight),sum(nominal_weight),sum(w_var)) # they are equal but u/d nominal are not equal
                    #
                    sumw, sumw2 = self.get_sumw_sumw2(df,w_var,y)
                    self.histos[f'{process}_{y}_{self.name}{ud_str}'] = {'sumw':sumw,
                                                                         'sumw2':sumw2}
                    
                
    def makeMCStatHist(self):
        for process in self.ids:
            for y in self.years:
                nom_hist = self.histos[f'{process}_{y}']
                stat_err = np.sqrt(nom_hist['sumw2'])
                mcstat_up, mcstat_down = nom_hist['sumw']+stat_err, nom_hist['sumw']-stat_err
                self.histos[f'{process}_{y}_{self.name}Up'] =   {'sumw' :mcstat_up,
                                                                 'sumw2':mcstat_up} # dont really care about sumw2 here
                self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw' :mcstat_down,
                                                                 'sumw2':mcstat_down} # dont really care about sumw2 here
                
    def handleQC(self):
        for process in self.ids:
            for y in self.years:
                nom     = np.float64(self.histos[f'{process}_{y}']['sumw'])
                nom_err = np.sqrt(np.float64(self.histos[f'{process}_{y}']['sumw2']))
                #print(f'{process}_{y}_{self.name}Up')
                up   = np.float64(self.histos[f'{process}_{y}_{self.name}Up']['sumw'])
                down = np.float64(self.histos[f'{process}_{y}_{self.name}Down']['sumw'])
                #UpRatio, DownRatio = np.divide(up_hist['sumw'],nom_hist['sumw']), np.divide(down_hist['sumw'],nom_hist['sumw'])
                #nom_sumw = np.where(nom_hist['sumw']==0, 0.00001 ,nom_hist['sumw'])
                # Step 1: kill_sys for bins where the stat err is larger than the nominal
                up   = np.where(nom<nom_err, nom, up)   
                down = np.where(nom<nom_err, nom, down)
                #UpRatio   = np.where(nom_sumw<np.sqrt(nom_hist['sumw2']), 1, up_hist['sumw']  / nom_sumw)
                #DownRatio = np.where(nom_sumw<np.sqrt(nom_hist['sumw2']), 1, down_hist['sumw']/nom_sumw)
                # Step 2: bad_sys
                #log_r_diff = abs(np.log(UpRatio)) - abs(np.log(DownRatio))
                #bad_sys = abs(log_r_diff) > 0.35
                #UpRatio, DownRatio = np.where(  bad_sys,np.where(log_r_diff > 0, 1/DownRatio, UpRatio), UpRatio), np.where(bad_sys,np.where(log_r_diff < 0, 1/UpRatio, DownRatio), DownRatio)
                # Step 3: one_sided_sys
                #one_sided_sys = (((UpRatio > 1) & (DownRatio > 1)) | ((UpRatio < 1) & (DownRatio < 1)))
                is_onesided = (( (up > nom) & (down > nom) ) | ( (up < nom) & (down < nom) ))
                #geo_mean = np.sqrt( np.divide( up, nom ) * np.divide( down, nom ) )
                geo_mean = np.sqrt( up * down )
                up , down = np.where(is_onesided, np.divide(up*nom, geo_mean) , up), np.where(is_onesided, np.divide(down*nom, geo_mean) , down) 
                #
                up, down = np.where((np.isinf(up)) | (np.isnan(up)), nom, up), np.where((np.isinf(down)) | (np.isnan(down)), nom, down) # handle cases where div by zero goes to inf
                
                # save hists 
                #self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = np.nan_to_num(UpRatio*nom_sumw)
                #self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = np.nan_to_num(DownRatio*nom_sumw)
                self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = up     # handle nan later
                self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = down   # handle nan later in TH1
        
    def handle_extraQC(self):
        for process in self.ids:
            for y in self.years:
                # get variation from merged bins for certain systematic
                # and distribute variation to finer bins, 
                # merged bins should be somewhat correlated in nature
                nom     = np.float64(self.histos[f'{process}_{y}']['sumw'])
                nom_err = np.sqrt(np.float64(self.histos[f'{process}_{y}']['sumw2']))
                #print(f'{process}_{y}_{self.name}Up')
                if re.search(r'(jms)|(jmr)', self.name) is not None and f'{process}_{y}_{self.name}Up' not in self.histos:
                    tmp_name = re.search(r'(jms)|(jmr)', self.name).group()+f'_{y}'
                    self.histos[f'{process}_{y}_{self.name}Up'] = {}     # need to initiate 
                    self.histos[f'{process}_{y}_{self.name}Down'] = {}
                else:
                    tmp_name = self.name
                up   = np.float64(self.histos[f'{process}_{y}_{tmp_name}Up']['sumw'])
                down = np.float64(self.histos[f'{process}_{y}_{tmp_name}Down']['sumw'])
                #
                # shape (3, 4, 4) pt,nn,sdm

                #if 'JES' in self.name or 'JEC' in self.name: # sum across NN per pt,sdM
                if 'hdamp' not in self.name and 'UE' not in self.name:
                    for i in range(nom.shape[-1]): # for sdM
                        for j in range(nom.shape[0]): # for pt
                            if 'fsr' in self.name : # merge across NN bins
                                up[j,:,i]   = nom[j,:,i] * np.nansum(up[j,:,i])/np.nansum(nom[j,:,i])
                                down[j,:,i] = nom[j,:,i] * np.nansum(down[j,:,i])/np.nansum(nom[j,:,i])
                            else : # dont merge variation in  first 2 bins
                                up[j,2:,i]   = nom[j,2:,i] * np.nansum(up[j,2:,i])/np.nansum(nom[j,2:,i])
                                down[j,2:,i] = nom[j,2:,i] * np.nansum(down[j,2:,i])/np.nansum(nom[j,2:,i])


                else: # dedicated sample sys, # sum across sdm and NN
                    #for i in range(nom.shape[-1]):
                    for j in range(nom.shape[0]): # pt
                        up[j,:,:]   = nom[j,:,:] * np.nansum(up[j,:,:])/np.nansum(nom[j,:,:])
                        down[j,:,:]   = nom[j,:,:] * np.nansum(down[j,:,:])/np.nansum(nom[j,:,:])
                        # dont merge variation in first 2 bins
                        #up[j,2:,:]   = nom[j,2:,:] * np.nansum(up[j,2:,:])/np.nansum(nom[j,2:,:])
                        #down[j,2:,:] = nom[j,2:,:] * np.nansum(down[j,2:,:])/np.nansum(nom[j,2:,:])

                self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = up     # handle nan later
                self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = down   # handle nan later in TH1
        
                
if __name__ == '__main__':
    #
    # initialize datacard making process
    hist2d = DataCardShapes(pt_bins,sdM_bins,n_NN_bins=10, isblind=isblind)# isblind false : cut out signal sensitive bins, isblind true: all bins# syntax np2d[year][ith_ptbin]
    MakeDataCard(isblind=isblind).makeDC()

