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
#
import uproot
#import uproot_methods
#from ROOT import TFile, TDirectory, TH1F
#import coffea.hist as hist
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt


np.random.seed(1)

'''
datacard bins defined here
use datacardshapes class to create nn bins
'''
#
isblind = False # this is still "blinded" if true, uses NN score  up to .7
#
pt_bins  = [0,200,300,450]
sdM_bins = [50,80,105,145,200]
#[50, 75, 85, 100, 110, 120, 135, 150, 200]
hist2d = DataCardShapes(pt_bins,sdM_bins,n_NN_bins=7, isblind=isblind) # syntax np2d[year][ith_ptbin]
#
# initializes np.histogram2d function with defined bins
def get_sumw_sumw2(df, weights, year):
    #sumw,  _    = self.hist3d([df['NN'].to_numpy(), df['Zh_M'].to_numpy(), self.ptclip(df['Zh_pt'].to_numpy())], # clip at 500 for large pt
    #                          weights=weights)
    #sumw2, _    = self.hist3d([df['NN'].to_numpy(), df['Zh_M'].to_numpy(), self.ptclip(df['Zh_pt'].to_numpy())], # clip at 500 for large pt
    # weights=np.power(weights,2))
    ret_NN_sdM = (lambda df: (df['NN'].to_numpy(), df['Zh_M'].to_numpy()))
    sumw = np.array([
        hist2d[year][i_bin]( *ret_NN_sdM(df[df['pt_bin']==i_bin]), 
                             weights=weights[df['pt_bin']==i_bin])[0] 
        for i_bin in range(1,len(pt_bins))
    ])
    sumw2 = np.array([
        hist2d[year][i_bin]( *ret_NN_sdM(df[df['pt_bin']==i_bin]), 
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
    #pt_bins = [0,200,300,450] # [200,300,400]
    #print(pt_bins[1:])
    #pt_bin_dict = {'200':'lopt', '300':'hipt'}
    dc_dir = 'Higgs-Combine-Tool'
    #output_rootfile = 'input_'+str(pt_bins[-1])+'inc.root'
    tag = '' if len(sys.argv) < 2 else sys.argv[1]+'_'
    #sig and bkg variables
    weights = ['weight','genWeight','Stop0l_topptWeight', #'SAT_HEMVetoWeight_drLeptonCleaned',
               'lep_trig_eff_tight_pt', 
               'lep_sf',
               'BTagWeightLight','BTagWeightHeavy','puWeight']#,'PrefireWeight']
    weight_sys = [#'Stop0l_topptWeight_Up' ,'Stop0l_topptWeight_Down', # only for top_powheg samples
        'BTagWeightLight_Up', 'BTagWeightLight_Down',
        'BTagWeightHeavy_Up', 'BTagWeightHeavy_Down',
        'puWeight_Up','puWeight_Down',
        'pdfWeight_Up','pdfWeight_Down',
        'ISR_Up', 'ISR_Down',
        'FSR_Up','FSR_Down',
        'mu_r_Up','mu_r_Down',
        'mu_f_Up','mu_f_Down',
        'mu_rf_Up','mu_rf_Down',
        #'PrefireWeight_Up','PrefireWeight_Down',
        'lep_trig_eff_tight_pt_up','lep_trig_eff_tight_pt_down',
        'lep_sf_up','lep_sf_down']
        
    bkg_v = weights + weight_sys + ['NN','Zh_M','Zh_pt','process']
    sig_v = bkg_v + ['genZHpt']
    #
    accepted_sig  = [f'{s}{i}' for i in range(len(pt_bins)) for s in ['ttH','ttZ']]#['ttHbb','ttZbb']] # change new and not new ttzbb here
    #accepted_bkg  = ['ttX','TTBar','old_tt_bb','Vjets','other']
    accepted_bkg  = ['ttX','TTBar','tt_bb','tt_2b','Vjets','other']
    accepted_data = ['data_obs']
    
    
    def __init__(self, 
                 sig = cfg.Sig_MC+cfg.sig_sys_samples, 
                 bkg = cfg.Bkg_MC+cfg.bkg_sys_samples,  # does not include QCD
                 years = cfg.Years, 
                 isblind=True):
        # start getting signal and bkg 
        self.sig = sig
        self.bkg = bkg
        self.data = cfg.Data_samples
        self.years =  years
        #self.years = ['2018']
        self.isblind = isblind
        self.dc_bins = len(pt_bins[1:])
        #
        #self.getdata()
        self.getdatav2()
        #print(self.data_dict.keys())
        self.initialize_hists() # switch to np.histogramdd , also, might be nice to get systematic shape plots here
        self.initialize_roofile() # format in file : $CHANNEL_$PROCESS, $CHANNEL_$PROCESS_$SYSTEMATIC
        self.initialize_datacard()
        # add systematics to histos
        self.add_Systematics()
        #print(self.histos.keys())
        if self.isblind:
            self.data_to_pdata()
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
        pool = Pool()
        print('Starting signal worker')
        all_samples = self.sig+self.bkg+self.data
        results = pool.map(self.worker, all_samples) # to run this in parallel 
        pool.close()
        pool.join()
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
        self.process_sig(self.data_dict) # add new_ttZbb meging , need to rewrite a bit
        #self.process_bkg(data_dict,y) # for tt_bb merging 
            
    def worker(self, process):
        y = re.search(r'201\d',process).group()
        p_vars = None # i think i have to do this to flush the memory
        p_vars = self.sig_v if process in self.sig else self.bkg_v
        p_vars = p_vars + cfg.ana_vars[f'sysvars_{y}']
        if '2018' in process : p_vars = p_vars + ['SAT_HEMVetoWeight_drLeptonCleaned']
        if 'pow' in process  : p_vars = p_vars + ['Stop0l_topptWeight_Up' ,'Stop0l_topptWeight_Down','tt_type']
        if 'Data' in process : p_vars = ['NN','Zh_M','Zh_pt','process']
        # get important info for signal
        return self.updatedict(process.replace(f'_{y}', ''), p_vars, y)
        
    #
        
    def updatedict(self, p, v, y=''):
        sub_f_dir = 'data_files' if 'Data' in p else 'mc_files'
        if not os.path.exists(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl'): return 
        df = pd.read_pickle(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl').filter(items=v)
        if 'TTbb' in p : df = pd.concat([df,pd.read_pickle(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl').filter(regex=r'\w*_weight')],axis='columns') # to get special ttbb normalizations
        df['Zh_pt'].clip(pt_bins[0]+1,pt_bins[-1]+1, inplace=True)
        df['pt_bin'] = pd.cut(df['Zh_pt'], bins=pt_bins+[500],
                              labels=[i_bin for i_bin in range(len(pt_bins))])
        group = df[df['NN'] >= 0.0].groupby(by='process')
        # extract sys type (if any)
        sys = '' 
        if p in cfg.all_sys_samples: 
            sys =  '_'+p.split('_')[-1] # format [process_name]_[systype]
            if 'hdamp' in sys and 'TTbb' in p: sys = sys.replace('hdamp','hdamp_ttbb') 
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
            s_pre = re.search(r'\w*[H,Z,Zbb]',sig_name).group()
            #print(sig_name)
            if sig_name not in data_dict: continue
            df = data_dict[sig_name]
            for i,_ in enumerate(pt_bins[:-1]):
                # insert 1 after sig_goup name # format should be [siggroup]_[sys]_[year]
                new_sig_name = f"{s_pre}{i}{sig_name.replace(s_pre, '')}"
                data_dict[new_sig_name] = df[(df['genZHpt'] >= pt_bins[i]) & 
                                             (df['genZHpt'] < pt_bins[i+1])]
                #
            #data_dict[f'{s}{len(pt_bins) -1 }_{y}'] = df[df['genZHpt'] >= pt_bins[-1]]
            data_dict[f"{s_pre}{len(pt_bins) -1}{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= pt_bins[-1]]
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
            y = re.search(r'201\d',s).group() # better than previous method
            w = getZhbbWeight(v,y) if 'data' not in s else np.ones_like(v['NN'].to_numpy()) # which should be 1 to handle data
            sumw, sumw2 = get_sumw_sumw2(v, w, y)
            self.histos[s] = {'sumw':sumw, # dim 1,2,3 is zh_pt, nn, zh_m
                              'sumw2':sumw2}
        

    @t2Run
    def add_Systematics(self):
        # SIG = ttZbb[0,1,2,3] ttHbb[0,1,2,3]
        # BKG = ttX, TTBar, old_tt_bb, Vjets, other
        process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        all_mc = self.accepted_sig + self.accepted_bkg
        all_but_ttbb = self.accepted_sig + ['TTBar','ttX','Vjets','other']
        tth_sig = re.findall(r'ttH\d', ' '.join(self.accepted_sig))
        ttz_sig = re.findall(r'\w*ttZ\d', ' '.join(self.accepted_sig))
        ttbar_mc = ['TTBar','tt_bb', 'tt_2b']
        #ttbar_mc = ['TTBar','old_tt_bb']
        #new_ttz_sig = ['new_'+z for z in ttz_sig]
        #
        Systematic.set_dc_processes(self.dc_dict, process_line)
        self.write2dc(f'# Rate uncertainties\n')
        Systematic('lumi_2016', 'lnN',  all_mc, 1.025)
        Systematic('lumi_2017', 'lnN',  all_mc, 1.023)
        Systematic('lumi_2018', 'lnN',  all_mc, 1.025)
        # probably do xsec theo rates here
        # first signal
        Systematic('tth_ggpdf', 'lnN', tth_sig, 1.036)
        Systematic('ttz_ggpdf', 'lnN', tth_sig, 1.035)
        Systematic('tth_qsc' ,  'lnN', tth_sig, 1.058,0.908)
        Systematic('ttz_qsc'  , 'lnN', ttz_sig, 1.081,0.907) 
        # now background
        Systematic('ggpdf', 'lnN', ttbar_mc, 1.042)
        Systematic('qqpdf', 'lnN', ['ttX','Vjets','other'],
                   [p_norms.rate_unc['pdf']['ttX'],
                    p_norms.rate_unc['pdf']['Vjets'],
                    p_norms.rate_unc['pdf']['other']])# combine into one 
        #Systematic('qqpdf', 'lnN', ['Vjets'], 1.038)# combine into one
        #Systematic('qqpdf', 'lnN', ['other'], 1.050)# combine into one
        #Systematic('tt_qsc'   , 'lnN', ttbar_mc+['ttX'], [[1.024,0.965] for _ in ttbar_mc]+[1.300])
        Systematic('tt_qsc'   , 'lnN', ttbar_mc, [[1.024,0.965] for _ in ttbar_mc])
        Systematic('ttx_qsc'  ,  'lnN', ['ttX'], p_norms.rate_unc['QCD_scale']['ttX']) 
        Systematic('v_qsc'    , 'lnN', ['Vjets'],p_norms.rate_unc['QCD_scale']['Vjets'])#1.008, 0.996) # .821/1.24
        Systematic('other_qsc', 'lnN', ['other'],p_norms.rate_unc['QCD_scale']['other']) #1.05, 0.95)   # .898/1.12
        Systematic('tt2bxsec', 'lnN', ['tt_2b'], 1.5)
        #Systematic('tt2bxsec', 'lnN', ['old_tt_bb'], 1.5)
        # Shape Systatics
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Shape uncertainties \n')
        ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos)#, self.hist3d, self.ptclip)
        # when defining shape, must incluse whether it is a mcsta, scale, or up/down syst
        for y in self.years:
            #self.histos = ShapeSystematic(f'btg_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_Up','BTagWeight_Down').get_shape()
            self.histos = ShapeSystematic(f'btglf_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeightLight_Up','BTagWeightLight_Down').get_shape()
            self.histos = ShapeSystematic(f'btghf_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeightHeavy_Up','BTagWeightHeavy_Down').get_shape()
            self.histos = ShapeSystematic(f'lepsf_{y}', 'shape', 'up/down', all_mc, 1, 'lep_sf_up','lep_sf_down').get_shape()
            self.histos = ShapeSystematic(f'trigeff_{y}', 'shape', 'up/down', all_mc, 1, 'lep_trig_eff_tight_pt_up','lep_trig_eff_tight_pt_down').get_shape()
        #
        self.histos = ShapeSystematic(f'pref_2016', 'shape', 'up/down', all_mc, 1, 'PrefireWeight_Up' ,'PrefireWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pref_2017', 'shape', 'up/down', all_mc, 1, 'PrefireWeight_Up' ,'PrefireWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'toppt', 'shape', 'up/down', ttbar_mc, 1, 'Stop0l_topptWeight_Up' ,'Stop0l_topptWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'isr', 'shape', 'ps', ['TTBar'], 1, 'ISR_Up','ISR_Down').get_shape()
        self.histos = ShapeSystematic(f'fsr', 'shape', 'ps', ['TTBar'], 1, 'FSR_Up','FSR_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_r', 'shape', 'normup/down', all_but_ttbb, 1, 'mu_r_Up','mu_r_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_f', 'shape', 'normup/down', all_but_ttbb, 1, 'mu_f_Up','mu_f_Down').get_shape()
        self.histos = ShapeSystematic(f'isr_ttbb', 'shape', 'ps', ['tt_bb', 'tt_2b'], 1, 'ISR_Up','ISR_Down').get_shape()
        self.histos = ShapeSystematic(f'fsr_ttbb', 'shape', 'ps', ['tt_bb', 'tt_2b'], 1, 'FSR_Up','FSR_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_r_ttbb', 'shape', 'normup/down', ['tt_bb', 'tt_2b'], 1, 'mu_r_Up','mu_r_Down').get_shape()
        self.histos = ShapeSystematic(f'mu_f_ttbb', 'shape', 'normup/down', ['tt_bb', 'tt_2b'], 1, 'mu_f_Up','mu_f_Down').get_shape()
        # redundant #self.histos = ShapeSystematic(f'mu_rf', 'shape', 'normup/down', all_mc, 1, 'mu_rf_Up','mu_rf_Down').get_shape()
        #
        #self.histos = ShapeSystematic(f'pdf_ttz', 'shape', 'up/down', ttz_sig, 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pdf', 'shape', 'normup/down', all_but_ttbb, 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pdf_ttbb', 'shape', 'normup/down', ['tt_bb','tt_2b'],                    1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        self.histos = ShapeSystematic(f'pu', 'shape', 'up/down', all_mc, 1, 'puWeight_Up','puWeight_Down').get_shape()
        # These shapes are already computed, just need to add to datacard
        self.histos = ShapeSystematic('ak4JER', 'shape', 'qconly', all_mc, 1).get_shape()
        self.histos = ShapeSystematic('ak4JES', 'shape', 'qconly', all_mc, 1).get_shape()
        self.histos = ShapeSystematic('ak8JER', 'shape', 'qconly', all_mc, 1).get_shape()
        self.histos = ShapeSystematic('ak8JES', 'shape', 'qconly', all_mc, 1).get_shape() 
        self.histos = ShapeSystematic('UE',     'shape', 'qconly', ['TTBar'], 1).get_shape()
        #self.histos = ShapeSystematic('erdOn', 'shape',  ['TTBar'], 1) # not working with just one shape at the moment
        self.histos = ShapeSystematic('hdamp', 'shape', 'qconly', ['TTBar'], 1).get_shape()
        self.histos = ShapeSystematic('hdamp_ttbb', 'shape', 'qconly', ['tt_bb', 'tt_2b'], 1).get_shape()
        #
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Float tt_bb normalization\n') 
        #self.write2dc('CMS_ttbbnorm rateParam * tt_*b 1 [-5,5]\n')
        Systematic('CMS_ttbbnorm', 'lnN', ['tt_bb','tt_2b'], 10)
        self.write2dc(100*'-'+'\n')
        self.write2dc('# MC Stats uncertainties\n') 
        #self.histos = ShapeSystematic(f'mcstat','shape','mcstat', all_mc, 1).get_shape()
        self.write2dc('* autoMCStats 10 0  1\n') 
        #

        

    def write2dc(self, str2write):
        for y in self.years:
            self.dc_dict[y].write(str2write)
            
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
                temp_dict = {'sumw' : v['sumw'][pt_bin,:,:].flatten()}#* (1 if y != '2017' else 3.3032)}#2.2967)} # to just scale to full run2
                hist_name = f'Zhpt{pt_bin+1}_{process}{sys}'
                if 'sumw2' in v:
                    temp_dict['sumw2'] = v['sumw2'][pt_bin,:,:].flatten()
                self.roo_dict[y][hist_name] = export1d(temp_dict, hist_name)

    def initialize_datacard(self):
        # creat 3 data cards , 1 per year
        dc_dict = {y: open(f'datacard_{self.tag}{y}.txt', 'w') for y in self.years}
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
            roo_name = f'datacard_{self.tag}{y}.root'
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
        self.years     = re.findall(r'201\d', name)
        if len(self.years) == 0: self.years = cfg.Years
        self.value    = value if type(value) is not list else {i:v for i,v in zip(self.ids,value)}
        self.optvalue = optvalue
        self.info     = '' if info is None else info
        #
        if self.datacard is not None: 
            for year in self.years:
                self.datacard[year].write(self.get_DC_line()) # write to datacard file upon instance creation

    @property
    def line(self):
        return '{0:14} {1:6}'.format(self.name,self.stype)

    @classmethod
    def set_dc_processes(cls,datacard, processes):
        cls.datacard = datacard
        cls.allowed_processes = processes

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

    def __init__(self, name, stype, subtype, process_ids, value, up=None, down=None, info=None):
        super().__init__(name, stype, process_ids, value, info)
        self.up     = up
        self.down   = down
        self.subtype = subtype

    def get_shape(self):
        fun_dict = {'mcstat' :self.makeMCStatHist,
                    #'scale'  :self.makeScaleHist,
                    'up/down'    :self.makeUpDownHist,
                    'normup/down':self.makeUpDownHist,
                    'ps'     :self.makeUpDownHist,#self.handlePSUpDown,
                    'erdOn'  :self.handleErdOn} # unique case where PS weights for 2017 tth and PS weight for 2016 and 2017 need to be handled 
        if self.subtype in fun_dict: fun_dict[self.subtype]()
        #
        #self.handleQC()
        #
        return self.histos

    @classmethod
    def set_df_histos_histfuncs(cls, data_dict, hist_dict):# hist_func, ptclip):
        cls.data_dict = data_dict
        cls.histos = hist_dict # to basically add histos dict to scope
        #cls.hist3d = hist_func
        #cls.ptclip = ptclip

    #def get_sumw_sumw2(self, df, weights):
            #sumw,  _    = self.hist3d([df['NN'].to_numpy(), df['Zh_M'].to_numpy(), self.ptclip(df['Zh_pt'].to_numpy())], # clip at 500 for large pt
            #                          weights=weights.to_numpy())
            #sumw2, _    = self.hist3d([df['NN'].to_numpy(), df['Zh_M'].to_numpy(), self.ptclip(df['Zh_pt'].to_numpy())], # clip at 500 for large pt
            #                          weights=np.power(weights.to_numpy(),2))
    #        return sumw, sumw2

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
                w_nom  = self.up.replace('_'+self.up.split('_')[-1],'') # format should be weightName_up/down
                #print(self.name, process)
                if w_nom == 'pdfWeight' or 'mu_' in w_nom or self.subtype == 'ps': 
                    df[w_nom] = 1.0
                weight =  getZhbbWeight(df,y)
                for ud_str, ud in zip(['Up','Down'],[self.up,self.down]):
                    nominal_weight = weight
                    #if 'new_tt_' in process and self.subtype is 'ps':
                        #nominal_weight = (weight/df['weight']) * df[f'{ud}_weight'] # to insert correct event weight
                    #    if 'mu_' in w_nom: df[w_nom] = df[f'{ud}_weight'] # to fix issue with tt_bb
                    w_var = nominal_weight*df[ud]/df[w_nom]
                    if 'norm' in self.subtype or ('tt_' in process and self.subtype == 'ps'): 
                        w_var = w_var * p_norms.shape_unc_norms[y][re.sub(r'\d$','', process)][ud]#sum(weight)/sum(w_var) # to keep the nominal normalization
                    #
                    #if 'new_tt_' in process and ('mu_' in w_nom or self.subtype is 'ps'):
                    #    print(process, y, ud, ud_str)
                    #    print(sum(weight),sum(nominal_weight),sum(w_var)) # they are equal but u/d nominal are not equal
                    #
                    sumw, sumw2 = get_sumw_sumw2(df,w_var,y)
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
                nom_hist  = self.histos[f'{process}_{y}']
                up_hist   = self.histos[f'{process}_{y}_{self.name}Up']
                down_hist = self.histos[f'{process}_{y}_{self.name}Down']
                nom_sumw = np.where(nom_hist['sumw']==0, 0.00001 ,nom_hist['sumw'])
                # Step 1: kill_sys
                UpRatio   = np.where(nom_sumw<np.sqrt(nom_hist['sumw2']), 1, up_hist['sumw']  / nom_sumw)
                DownRatio = np.where(nom_sumw<np.sqrt(nom_hist['sumw2']), 1, down_hist['sumw']/nom_sumw)
                print(UpRatio)
                #print(DownRatio)
                # Step 2: bad_sys
                log_r_diff = abs(np.log(UpRatio)) - abs(np.log(DownRatio))
                bad_sys = abs(log_r_diff) > 0.35
                UpRatio, DownRatio = np.where(  bad_sys,np.where(log_r_diff > 0, 1/DownRatio, UpRatio), UpRatio), np.where(bad_sys,np.where(log_r_diff < 0, 1/UpRatio, DownRatio), DownRatio)
                # Step 3: one_sided_sys
                one_sided_sys = (((UpRatio > 1) & (DownRatio > 1)) | ((UpRatio < 1) & (DownRatio < 1)))
                UpRatio, DownRatio = np.where(one_sided_sys, UpRatio/np.sqrt(UpRatio*DownRatio),UpRatio), np.where(one_sided_sys, DownRatio/np.sqrt(UpRatio*DownRatio),DownRatio)
                # save hists 
                self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = np.nan_to_num(UpRatio*nom_sumw)
                self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = np.nan_to_num(DownRatio*nom_sumw)
        

class DCParam :
    '''
    Handles different parameter types
    to be added to datacard
    '''    
    datacard     = None
    #
    def __init__(self, name, param_type, process, i_val=None, sigma=None, up=None, down=None, write=True):
        self.name       = name
        self.param_type = param_type
        self.process    = process
        self.i_val      = i_val
        self.sigma      = sigma
        self.up         = up
        self.down       = down
        self.param_dict = {'param'    :self.get_param_line(),
                           'flatParam':self.get_param_line(),
                           'extArg'   :self.get_extarg_line(),
                           'ratParam' :self.get_rateparam_line()
        }
        if write:  self.write_to_dc()
    
    @classmethod
    def set_current_dc(cls, datacard):
        cls.datacard = datacard
    
    def get_param_line(self):
        _line = f'{self.name:8} {self.param_type:14} {self.i_val:14} {self.sigma:14}\n'
        return _line
    def get_extarg_line(self):
        _line = f'{self.name:8} {self.param_type:14} {self.i_val:6} [{self.up},{self.down}]\n'
        return _line
    def get_rateparam_line(self):
        _line = f'{self.name:8} {self.param_type:14} *\t {self.process:14} {self.i_val:6} [{self.up},{self.down}]\n'
        return _line
    def write_to_dc(self):
        self.datacard.write(self.param_dict[self.param_type])

class EFTParam(DCParam):
    '''
    Inheiriting DCParam class to handle the task
    of adding relevent EFT parameter rates
    including WC args, P Q R rates and 
    overall EFT scale factor per gen bin
    
    Must specify ttH or ttZ
    '''
    
    def __init__(self, sig_type):
        self.sig_type = sig_type
        self.pqr_df   = pd.from_pickle('pqr_df.pkl')
        self.wc_range = pd.from_pickle('wc_minmax_df.pkl')
        self.datacard = super().datacard
        # add WC ext args using DCParam class
        # then take care of ther remaining 55 parameters, 1 SM, 9 Q, 45 P

    def add_WC_to_DC(self):
        #dcparam = super().__init__(args go here for WC, get WC from df)
        pass
        
        
                
if __name__ == '__main__':
    #
    # initialize datacard making process
    MakeDataCard(isblind=isblind)    

