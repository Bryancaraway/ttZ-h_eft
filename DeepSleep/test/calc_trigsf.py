##### CALC LEPTON TRIG EFF SF  #####
### Written by: Bryan Caraway    ###
####################################
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import uproot
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
#import functools
import config.ana_cff as cfg
from config.sample_cff import sample_cfg, process_cfg
from modules.AnaDict import AnaDict
from modules.getdata import getData
from lib.fun_library import clop_pear_ci


class Calc_LepEffSF :
    '''
    calc trig eff sf per 
    single muon and single electron
    channel per year using the opposite
    flavor as a reference trigger
    '''
    of_lep = {
        'Electron':'Muon',
        'Muon':'Electron',
    }
    def __init__(self,channel, year):
        self.channel = channel # Electron or Muon
        self.year = year
        # load data and mc that passes reference trigger
        self.data_df = self.getData('Data_Single'+self.of_lep[self.channel])
        self.mc_df   = self.getData('TTTo2L2Nu')
        # need to load and add sf to mc_df
        self.add_refTrigSF()
        

    def add_refTrigSF(self,period):
        filesf_dict = {
            'Electron':{'2016':(lambda : json.load(open(".json"))),
                        '2017':(lambda : json.load(open(".json"))),
                        '2018':(lambda : json.load(open(".json"))),
                    },
            'Muon':    {'2016_AtoF'  : (lambda : uproot.open("EfficienciesAndSF_RunBtoF.root")['IsoMu24_OR_IsoTkMu24_PtEtaBins']['pt_abseta_ratio']),
                        '2016_GtoH'  : (lambda : uproot.open("EfficienciesAndSF_Period4.root")['IsoMu24_OR_IsoTkMu24_PtEtaBins']['pt_abseta_ratio']),
                        '2017'       : (lambda : uproot.open("EfficienciesAndSF_RunBtoF_Nov17Nov2017.root")['IsoMu27_PtEtaBins']['pt_abseta_ratio']),
                        '2018_before': (lambda : uproot.open("EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root")['IsoMu24_PtEtaBins']['pt_abseta_ratio']),
                        '2018_after' : (lambda : uproot.open("EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root")['IsoMu24_PtEtaBins']['pt_abseta_ratio']),
                    }
            }
        ref_trigger = self.filesf_dict[self.of_lep[self.channel]][period]()
        #### i am here ####
            

    def getData(self,sample):
        input_file = cfg.postSkim_dir+f"{self.year}/Trig_{sample_cfg[sample]['out_name']}/{sample}.pkl"
        isData   = 'Data' in sample
        gD_out = AnaDict.read_pickle(f'{input_file}')
        return gD_out['events'][self.pass_refTrig(gD_out['events'])]
        #if isData:
        #    self.data_df = pd.concat([self.data_df,gD_out['events']], axis='rows', ignore_index=True)
    def pass_refTrig(self,df):
        ref_trig_dict = {
            'Electron':{'2016': (lambda x : ((x['HLT_IsoMu24']) | 
                                             (x['HLT_IsoTkMu24']))),
                        '2017': (lambda x : (x['HLT_IsoMu27'])),
                        '2018': (lambda x : (x['HLT_IsoMu24'])),
                    },
            'Muon'    : cfg.hlt_path['electron'],
        }
        pass_reftrig = ref_trig_dict[self.channel][self.year](df)
        return pass_reftrig

if __name__ == '__main__':
    Calc_LepEffSF('Electron','2016')
