########################
### Skim data        ###
### for trigger      ###
### efficiecny sf    ###
### calc             ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#
import os
import sys
import json
import numpy as np
import pandas as pd
#from modules.metaSkim import SkimMeta
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import config.ana_cff as cfg
from modules.Skim import Skim


class TrigSkim(Skim) :
    '''
    This class borrows the functions
    and obj selection definitions from
    the super class: SKim
    --> with changes to event selection
    reqs, and lepton object selection
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)

    def get_event_selection(self): # after objects have been defined
        return ( (self.jets['Jet_pt'].counts >= 4) &
                 (self.events['MET_pt'] > 20) &
                 (self.electrons['Electron_pt'].counts == 1) & 
                 (self.muons['Muon_pt'].counts == 1) &
                 self.get_MET_filter() &
                 (self.pass_goldenjson() == 1 if self.isData else True) &
                 (
                     np.where(self.events['run'] >= 319077, self.get_HEM_veto(), True) if 
                     (self.year == '2018' and self.isData) else True
                 )
             )
    def handle_lep_info(self): # assumes event filter already applied (1 lepton per event)
        for el_mu, elmu_df in zip(['Electron','Muon'],[self.electrons, self.muons]):
            self.events.update({
                f'{el_mu}_pt' : elmu_df[f'{el_mu}_pt'].flatten(),
                f'{el_mu}_eta': elmu_df[f'{el_mu}_eta'].flatten(),
                f'{el_mu}_phi': elmu_df[f'{el_mu}_phi'].flatten(),
                f'{el_mu}_mass': elmu_df[f'{el_mu}_mass'].flatten(),
                f'{el_mu}_sip3d': elmu_df[f'{el_mu}_sip3d'].flatten(),
                f'{el_mu}_charge': elmu_df[f'{el_mu}_charge'].flatten(),
            })

    def get_b_tag_eventmask(self):
        return (self.events['nBottoms'] >= 0)
        
    def get_skim(self):
        __out_dict = {}
        events = {**self.events,**self.hlt}
        __out_dict['events'] = pd.DataFrame.from_dict({k:v for k,v in events.items()})
        __out_dict['soft_e'] = self.soft_electrons
        __out_dict['soft_mu'] = self.soft_muons
        if not self.isData:
            __out_dict['gen'] = self.geninfo
            __out_dict['metaData'] = self.Meta.get_metadata()
        # close root file
        self.f.close()
        #
        return __out_dict

    def is_a_electron(self):
        return  cfg.lep_sel['nosip3d_electron'][self.year](self.electrons)

    def is_a_muon(self):
        return cfg.lep_sel['nosip3d_muon'](self.muons)

if __name__ == '__main__':
    year = '2016'
    golden_json=json.load(open(cfg.goodLumis_file[year]))
    test_file   = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2016/TTTo2L2Nu_2016/8D57A168-9993-A34B-A743-E04C47F32D1D_Skim_7.root'
    _ = TrigSkim(test_file, 'TTTo2L2Nu', year, isData=False, jec_sys=None, golden_json=golden_json)
    print(_.get_skim())
