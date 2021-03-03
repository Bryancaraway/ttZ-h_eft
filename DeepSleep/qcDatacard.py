#!/usr/bin/env python

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
from makeDatacard import MakeDataCard, Systematic, ShapeSystematic
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
target = ''
if len(sys.argv) > 1:
    target = sys.argv[1]
doNNcuts = False
if len(sys.argv) > 2: # means do nn cuts
    doNNcuts = True

#nn = 'NN'
nn = cfg.nn
    
# only targets listed here may be used
tbins_map = {
    #'Zh_M' :[50,80,105,145,200],
    #'Zh_pt':[200,300,450,np.inf],
    #'nBottoms_drLeptonCleaned':[1.5,2.5,3.5,4.5,10],

    # NN variables onward
    'outZH_b1_pt':[30,60,90,150,210,300,500,np.inf],
    'outZH_b2_pt':[30,60,90,150,210,np.inf],
    'outZH_b1_score':[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],
    'outZH_b2_score':[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],    
    'outZh_q1_pt':[30,60,90,150,210,300,500,np.inf],
    'outZh_q2_pt':[30,60,90,150,210,300,500,np.inf],
    'outZh_q1_btag':[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],
    'outZh_q2_btag':[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],    
    'outZH_b1_q_mindr':[0.4, 0.8, 1.2, 2.4, np.inf],
    'outZH_b2_q_mindr':[0.4, 0.8, 1.2, 2.4, np.inf],
    'outZH_q_q_dr_nearb1': [0.4, 0.8, 1.2, 2.4, 3.6, np.inf],
    'outZH_q_q_dr_nearb2': [0.4, 0.8, 1.2, 2.4, 3.6, np.inf],
    'outZH_qq_M_nearb1':[0,50,100,150,200,np.inf],
    'outZH_qq_M_nearb2':[0,50,100,150,200,np.inf],
    'outZH_b1q_M':[0,50,100,150,200,np.inf],
    'outZH_b1_qq_dr': [0.0, 1.2, 2.4, np.inf], 
    'outZH_b2_qq_dr': [0.0, 1.2, 2.4, np.inf], 
    'outZH_b1qq_M': [0,50,100,150,200,np.inf],
    'outZH_b2qq_M': [0,50,100,150,200,np.inf],
    'ZH_b1qq_dr':  [0.0, 0.8, 1.2, 2.4, 3.6, np.inf], 
    'ZH_b2qq_dr':  [0.0, 0.8, 1.2, 2.4, 3.6, np.inf], 
    'ZH_lbb1qq_dr': [0.0, 0.8, 1.2, 2.4, 3.6, np.inf], 
    'ZH_lbb2qq_dr': [0.0, 0.8, 1.2, 2.4, 3.6, np.inf], 
    'l_b2_mtb': [0,50,100,150,200,np.inf],
    'Zh_closeb_invM': [0,100,150,200,300,np.inf],  
    'n_ak8jets':[-0.5,0.5,1.5,2.5,np.inf],    
    'n_ak4jets':[3.5,4.5,5.5,6.5,7.5,np.inf], # NN
    'n_ak8_Zhbb':[-0.5,0.5,1.5,np.inf],
    'outZh_max_ak8sdM' : [50,100,150,200,np.inf],
    'outZh_b12_m' : [0,50,100,150,200,np.inf],
    'outZh_b12_dr': [0.4, 1.2, 2.4, 3.6, np.inf],
    'ht_b': [0,200,300,500,700,np.inf],
    'ht_outZh': [0,300,500,700,900,np.inf],
    'n_Zh_btag_sj':[-0.5,0.5,1.5,np.inf],
    'Zh_bestb_sj':[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],
    'Zh_worstb_sj':[0,.25,.5,.75,1.0],
    'nonZhbb_b1_dr':[0.8, 1.2, 1.6, 2.0, 2.4, np.inf],
    'nonZhbb_q1_dr':[0.8, 1.2, 1.6, 2.0, 2.4, np.inf],
    'inZhb_outZhb_dr': [0.0, 0.8, 1.2, 1.6, 2.0, 2.4, np.inf], 
    'Zh_l_dr':[0.8, 1.2, 1.6, 2.0, 2.4, np.inf],
    'Zh_l_invM_sd':[50,100,125,150,175,200,np.inf],
    'l_b1_invM':[50,100,125,150,175,200,np.inf],
    'l_b2_invM':[50,100,125,150,175,200,np.inf],
    'l_b1_dr':[0.4,0.8, 1.2, 1.6, 2.0, 2.4, np.inf],
    'l_b2_dr':[0.4, 1.2, 1.6, 2.0, 2.4, np.inf],
    'spher':[0,.1,.2,.3,.4,.5,.6,.7,.8,1.0],
    'aplan':[0,.025,.05,.075,.1,.15,.20,.25,.30,np.inf],
    'n_b_inZh':[-0.5,0.5,1.5,np.inf], # NN
    'n_q_inZh':[-0.5,0.5,1.5,np.inf],
    'n_q_outZh':[-0.5,0.5,1.5,2.5,3.5,4.5,np.inf],# NN
    'n_b_outZh':[1.5,2.5,3.5,4.5,np.inf], 
    'Zh_bbvLscore':[.8,.85,.9,1.0], # NN
    # ---
    # Control Plots
    # ---
    

}
if target not in tbins_map and __name__ == '__main__':
    raise KeyError(f"{target} not in 'tbins_map'!!!\nOnly the following are available: {list(tbins_map.keys())}")
    
def get_sumw_sumw2(df, weights, year):
    #ret_target = (lambda df: (df['NN'].to_numpy(), df['Zh_M'].to_numpy()))
    #sumw = np.array([
    #    hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
    #                         weights=weights[df['pt_bin']==i_bin])[0] 
    #    for i_bin in range(1,len(pt_bins))
    #])
    #sumw2 = np.array([
    #    hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
    #                         weights=np.power(weights[df['pt_bin']==i_bin],2))[0]
    #    for i_bin in range(1,len(pt_bins))
    #])
    sumw  = np.histogram( df[target].to_numpy(),
                          bins = tbins_map[target],
                          weights=weights)[0]
    sumw2 = np.histogram( df[target].to_numpy(),
                          bins = tbins_map[target],
                          weights=np.power(weights,2))[0]

    return sumw, sumw2
#

class MakeQCDataCard(MakeDataCard):
    '''
    Class to handle the overhead of creating,
    and writing the datacard as well as 
    the root file with shapes for template fit
    '''
    #file_dir = './files/'
    #
    #pt_bins = [0,200,300,450] # [200,300,400] # for gen Zhpt 
    pt_bins = [0,450] # [200,300,400] # for gen Zhpt 
    accepted_sig = ['ttH','ttZ']
    #dc_dir = 'Higgs-Combine-Tool'
    tag = target+'_'+ ( 'NNcuts_' if doNNcuts else '' )
        
    def __init__(self, 
                 sig = cfg.Sig_MC+cfg.sig_sys_samples, 
                 bkg = cfg.Bkg_MC+cfg.bkg_sys_samples,  # does not include QCD
                 years = cfg.Years, 
                 isblind=False):
        # start getting signal and bkg 
        super().__init__(sig,bkg,years,isblind,get_sumw_sumw2)
        self.sig = sig
        self.bkg = bkg
        self.data = cfg.Data_samples
        self.years =  years
        #self.years = ['2018']
        self.isblind = isblind
        self.dc_bins = 1 #len(tbins_map[target])#len(pt_bins[1:])
        self.bkg_v  = super().weights + super().weight_sys + ['process','sample',target]
        self.sig_v  = self.bkg_v + ['genZHpt']
        self.data_v = ['process',target]
        #
        #self.getdata()

    def makeQCDC(self):
        super().getdatav2()
        #super().process_sig(self.data_dict)
        #print(self.data_dict.keys())
        super().initialize_hists() 
        super().initialize_roofile() # format in file : $CHANNEL_$PROCESS, $CHANNEL_$PROCESS_$SYSTEMATIC
        self.initialize_datacard()
        # add systematics to histos
        self.setup_Systematics()
        super().add_Systematics()
        #super().add_uncorrSystematics()
        #print(self.histos.keys())
        if self.isblind:
            self.data_to_pdata()
        self.fill_roofile()
        super().close_roofile()
        super().close_dc()
        
    def updatedict(self, p, v, y=''):
        if doNNcuts:
            v = v +[nn]
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

        #df['Zh_pt'].clip(pt_bins[0]+1,pt_bins[-1]+1, inplace=True)
        #df['pt_bin'] = pd.cut(df['Zh_pt'], bins=pt_bins+[500],
        #                      labels=[i_bin for i_bin in range(len(pt_bins))])
        if doNNcuts:
            group = df[df[nn] >= 0.0].groupby(by='process')
        else:
            group = df.groupby(by='process')
        if 'Data' not in p:
            for n,g in group:
                if n not in ['ttH','ttZ'] + self.accepted_bkg: 
                    df.drop(g.index, inplace=True)
        # extract sys type (if any)
        sys = '' 
        data_dict = {}
        if p in cfg.all_sys_samples: 
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

            
    def setup_Systematics(self):
        process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        Systematic.set_dc_processes(self.dc_dict, process_line)
        ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos, self.get_sumw_sumw2)
        

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
            #for ibin in range(v['sumw'].shape[0]):
            #    if pt_bin == -1:
            #        to_flat = (lambda a: self.merge_last_mbin(pt_bin,a))
            #    else:
            #        to_flat = (lambda a : a[pt_bin,:,:].flatten())
            #    temp_dict = {'sumw' : to_flat(v['sumw'])}#* (1 if y != '2017' else 3.3032)}#2.2967)} # to just scale to full run2
                #hist_name = f'Zhpt{pt_bin+1}_{process}{sys}'
            hist_name = f'{target}_{process}{sys}'
            temp_dict = {'sumw' :v['sumw'],
                         'sumw2':v['sumw2']}
            #    if 'sumw2' in v:
            #        temp_dict['sumw2'] = to_flat(v['sumw2'])
            self.roo_dict[y][hist_name] = export1d(temp_dict, hist_name)


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
                #f"{'bin':20}{' '.join(['ch'+str(i+1) for i in range(self.dc_bins)])}\n",
                f"{'bin':20}{' '.join([target for i in range(self.dc_bins)])}\n",
                f"{'observation':20}{self.dc_bins*'-1  '}\n",
                100*'-'+'\n',
                #f"{'bin':20}{' '.join(['ch'+str(i+1) for i in range(self.dc_bins) for _ in range(len(self.accepted_sig + self.accepted_bkg))])}\n",
                f"{'bin':20}{' '.join([target for i in range(self.dc_bins) for _ in range(len(self.accepted_sig + self.accepted_bkg))])}\n",
                f"{'process':20}{' '.join(s for _ in range(self.dc_bins) for s in self.accepted_sig+self.accepted_bkg)}\n",
                f"{'process':20}{' '.join(str(i) for _ in range(self.dc_bins) for i in range(-len(self.accepted_sig)+1, len(self.accepted_bkg)+1))}\n",
                f"{'rate':20}{' '.join(str(-1) for _ in range(self.dc_bins) for s in self.accepted_sig + self.accepted_bkg)}\n",
                100*'-'+'\n'])
        #
        self.dc_dict = dc_dict

                
if __name__ == '__main__':
    #
    # initialize datacard making process
    MakeQCDataCard(isblind=False).makeQCDC()

