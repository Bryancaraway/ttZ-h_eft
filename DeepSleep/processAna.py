########################
### process code     ###
### for analysis     ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#

import sys
import os
import math
import pickle
import operator
import re
import time
from collections import defaultdict
from numba import njit, prange
import concurrent.futures
#
import deepsleepcfg as cfg
#import kinematicFit as kFit
import fun_library as lib
from fun_library import fill1e, fillne, deltaR, deltaPhi, invM, calc_mtb
#
import numpy as np
np.random.seed(0)
import pandas as pd
from itertools import combinations
##

class processAna :
    '''
    callable object designed to process data that has been extracted from
    root files
    
    the output are .pkl files ready for histogram analysis and Combine Tool
    '''
    # Allowed class variables
    files    = None
    samples  = None
    isData   = False
    isSignal = False
    isttbar  = False
    #
    lc     = '_drLeptonCleaned'
    #
    gen_df = None
    jet_df = None 
    val_df = None
    rtc_df = None
    #
    ak4_dR = 1.4
    ak8_dR = 0.8

    def __init__(self, kwargs):
        [setattr(self,k,v) for k,v in kwargs.items() if k in processAna.__dict__.keys()]
        self.ak4_df = self.jet_df.loc['AK4_Variables'].dropna(how='all',axis='columns') 
        self.ak8_df = self.jet_df.loc['AK8_Variables'].dropna(how='all',axis='columns')     
        del self.jet_df
        self.process_data()

    def process_data(self):
        # to contain lepCleaned_v2, match_gen, and ZHbbAna
        #self.lepCleaned_v2()
        self.match_gen_tt()

    def lepCleaned_v2(self):
        # Clean ak4 and ak8 jets within dR of 0.4 and 0.8 of Lepton
        lep_eta = self.val_df['Lep_eta'].to_numpy()
        lep_phi = self.val_df['Lep_phi'].to_numpy()
        iter_   = zip([self.ak4_df,self.ak8_df],['Jet','FatJet'], [self.ak4_dR, self.ak8_dR])
        for j_df, j_name, eta_cut in iter_:
            j_eta, j_phi = j_df[j_name+'_eta'+self.lc].unstack(), j_df[j_name+'_phi'+self.lc].unstack()
            j_keys = re.findall(j_name+r'_\w+', ' '.join(j_df.columns.values))
            lep_j_dR = deltaR(lep_eta,lep_phi, j_eta, j_phi).stack()
            j_df.loc[:,j_keys] = j_df.loc[:,j_keys].where(lep_j_dR > eta_cut, np.nan, axis='columns')
        #
    def match_gen_tt(self):
        # match tt or ttZ/h gen particles to recontructed objects
        print(self.gen_df.unstack())
        gen_ids = self.gen_df['GenPart_pdgId'].unstack(fill_value=999).to_numpy()
        gen_mom = self.gen_df['GenPart_genPartIdxMother'].unstack(fill_value=-2).to_numpy()
        gen_mom_ids = np.take_along_axis(gen_ids, gen_mom, axis = 1)
        is_tt_bb = (np.sum((abs(gen_ids) == 5) & (abs(gen_mom_ids) > 6) & (gen_mom_ids != 999), axis=1) >= 2)
        self.val_df['tt_bb'] = is_tt_bb
        
    def match_gen_sig(self):
        # match tt or ttZ/h gen particles to recontructed objects
        pass

    def match_gen_lep(self):
        lep_eta = self.val_df['Lep_eta'].to_numpy()
        lep_phi = self.val_df['Lep_phi'].to_numpy()
        #
        
        pass

    def ZHbbAna(self):
        pass


if __name__ == '__main__':
    #
    # will need to open pkl files for testing
    print('Reading Files...')
    dir_ = 'test_files'
    jet_df = pd.read_pickle(dir_+'/TTZH_2017_jet.pkl')
    val_df = pd.read_pickle(dir_+'/TTZH_2017_val.pkl')
    gen_df = pd.read_pickle(dir_+'/TTZH_2017_gen.pkl')
    print('Processing data...')
    process_ana_dict = {'jet_df':jet_df , 'val_df':val_df, 'gen_df':gen_df}
    processAna(process_ana_dict)
