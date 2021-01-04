########################
### pre Skim         ###
### for analysis     ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#

import uproot
import os
import re
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import time
from collections import defaultdict
from numba import njit, prange
import concurrent.futures
import multiprocessing
executor = concurrent.futures.ThreadPoolExecutor()
import functools
from glob import glob
#
import config.ana_cff as cfg
from config.sample_cff import sample_cfg
from lib.fun_library import fillne, t2Run
from modules.AnaDict import AnaDict
from modules.AnaVars import AnaVars
#
import numpy as np
np.random.seed(0)
import pandas as pd


class PreSkim :
    '''
    only run on MC
    get 4 main items from sample
    1) event normalization  (all) 
    2) varied normalization (all)
    3) varied normalization w/resp tt+bb (only ttbar and ttbb)
    4) r for btagSFweight per np, and per jet multi
    Store this info in MetaData.json
    '''
    isData = False
    outDir = '' # add later */*postSkim
    nanoAODv7_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/'
    sample  = None
    jec_sys = None
    min_njets    = 4
    min_nfatjets = 1
    min_nbjets   = 2
    
    def __init__(self, files, sample, year, isttbar=False, isttbb=False):
        self.sample = sample
        self.year   = year
        self.isttbar = isttbar
        self.isttbb  = isttbb
        self.metadata = {'sample':sample, 'year':year, 'xs': sample_cfg[sample]['xs'], 'kf': sample_cfg[sample]['kf']}
        self.files    = files#glob(f'{self.nanoAODv7_dir}/{year}/{sample}_{year}/*.root')
        self.pool     = multiprocessing.Pool()
        #
        if self.isttbar or self.isttbb:
            self.ttbar_ttbb_normalization()
        else:
            self.event_var_normalization()

    ## main functions    
    
    def event_var_normalization(self): # 1,2,4
        result_list = np.array(self.pool.map(self.event_worker, self.files))
        results =  [sum(result_list[:,i]) for i in range(26)]
        self.metadata.update({
            'p_events': results[0], 'n_events': results[1],
            'tot_events': results[2], 
            'mur_up_tot':results[3],'mur_down_tot':results[4],
            'muf_up_tot':results[5],'muf_down_tot':results[6],
            'murf_up_tot':results[7],'murf_down_tot':results[8],
            'pdf_up_tot':results[9],'pdf_down_tot':results[10],
            #
            'r_nom': results[11]/results[2],
            'r_jes_up': results[12]/results[2],'r_jes_down': results[13]/results[2],
            'r_lf_up': results[14]/results[2],'r_lf_down': results[15]/results[2],
            'r_hf_up': results[16]/results[2],'r_hf_down': results[17]/results[2],
            'r_lfstats1_up': results[18]/results[2],'r_lfstats1_down': results[19]/results[2],
            'r_lfstats2_up': results[20]/results[2],'r_lfstats2_down': results[21]/results[2],
            'r_hfstats1_up': results[22]/results[2],'r_hfstats1_down': results[23]/results[2],
            'r_hfstats2_up': results[24]/results[2],'r_hfstats2_down': results[25]/results[2],

            
        })

    def ttbar_ttbb_normalization(self): # 3,4
        result_list = np.array(self.pool.map(self.ttbb_worker, self.files))
        results =  [sum(result_list[:,i]) for i in range(41)]
        self.metadata.update({
            'tot_events': results[0], 
            'mur_up_tot':results[1],'mur_down_tot':results[2],
            'muf_up_tot':results[3],'muf_down_tot':results[4],
            'murf_up_tot':results[5],'murf_down_tot':results[6],
            'isr_up_tot':results[7],'isr_down_tot':results[8],
            'fsr_up_tot':results[9],'fsr_down_tot':results[10],
            'pdf_up_tot':results[11],'pdf_down_tot':results[12],
            #
            'ttbb_events': results[13], 
            'mur_up_ttbb':results[14],'mur_down_ttbb':results[15],
            'muf_up_ttbb':results[16],'muf_down_ttbb':results[17],
            'murf_up_ttbb':results[18],'murf_down_ttbb':results[19],
            'isr_up_ttbb':results[20],'isr_down_ttbb':results[21],
            'fsr_up_ttbb':results[22],'fsr_down_ttbb':results[23],
            'pdf_up_ttbb':results[24],'pdf_down_ttbb':results[25],
            #
            'r_nom': results[26]/results[0],
            'r_jes_up': results[27]/results[0],'r_jes_down': results[28]/results[0],
            'r_lf_up': results[29]/results[0],'r_lf_down': results[30]/results[0],
            'r_hf_up': results[31]/results[0],'r_hf_down': results[32]/results[0],
            'r_lfstats1_up': results[33]/results[0],'r_lfstats1_down': results[34]/results[0],
            'r_lfstats2_up': results[35]/results[0],'r_lfstats2_down': results[36]/results[0],
            'r_hfstats1_up': results[37]/results[0],'r_hfstats1_down': results[38]/results[0],
            'r_hfstats2_up': results[39]/results[0],'r_hfstats2_down': results[40]/results[0],

        })

    def get_metadata(self):
        return self.metadata
    

    ## worker functions 
    @staticmethod
    def event_worker(roo):
        with uproot.open(roo) as f:
            t = f['Events']
            gw = t.array('genWeight')#, executor=executor, blocking=True)
            scale = t.array('LHEScaleWeight').pad(9).fillna(1) * np.sign(gw)
            pdf_up   = t.array('pdfWeight_Up') * np.sign(gw)
            pdf_down = t.array('pdfWeight_Down') * np.sign(gw)
            #
            sc_tot = [ sum( scale[:,i] ) for i in range(9)]
            mur_up_tot, mur_down_tot   = sc_tot[7], sc_tot[1]
            muf_up_tot, muf_down_tot   = sc_tot[5], sc_tot[3]
            murf_up_tot, murf_down_tot = sc_tot[8], sc_tot[0]
            #
            pdf_up_tot, pdf_down_tot = sum(pdf_up), sum(pdf_down)
            #
            p_count= sum(gw>=0.0)
            n_count= sum(gw<0.0)
            tot_count = p_count - n_count
            #
            bt_nom = sum(t.array('BTagWeight') * np.sign(gw))
            bt_jes_up, bt_jes_down = sum(t.array('BTagWeight_jes_up') * np.sign(gw)), sum(t.array('BTagWeight_jes_down') * np.sign(gw))
            bt_lf_up, bt_lf_down = sum(t.array('BTagWeight_lf_up') * np.sign(gw)),    sum(t.array('BTagWeight_lf_down') * np.sign(gw))
            bt_hf_up, bt_hf_down = sum(t.array('BTagWeight_hf_up') * np.sign(gw)),    sum(t.array('BTagWeight_hf_down') * np.sign(gw))
            bt_lfstats1_up, bt_lfstats1_down = sum(t.array('BTagWeight_lfstats1_up') * np.sign(gw)), sum(t.array('BTagWeight_lfstats1_down') * np.sign(gw))
            bt_lfstats2_up, bt_lfstats2_down = sum(t.array('BTagWeight_lfstats2_up') * np.sign(gw)), sum(t.array('BTagWeight_lfstats2_down') * np.sign(gw))
            bt_hfstats1_up, bt_hfstats1_down = sum(t.array('BTagWeight_hfstats1_up') * np.sign(gw)), sum(t.array('BTagWeight_hfstats1_down') * np.sign(gw))
            bt_hfstats2_up, bt_hfstats2_down = sum(t.array('BTagWeight_hfstats2_up') * np.sign(gw)), sum(t.array('BTagWeight_hfstats2_down') * np.sign(gw))
            #
            return [p_count, n_count, # 0-1
                    tot_count,        # 2
                    mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot, # 3-8
                    pdf_up_tot, pdf_down_tot, # 9-10
                    bt_nom, # 11
                    bt_jes_up, bt_jes_down, bt_lf_up, bt_lf_down, bt_hf_up, bt_hf_down, # 12-17
                    bt_lfstats1_up, bt_lfstats1_down, bt_lfstats2_up, bt_lfstats2_down,  # 18-21
                    bt_hfstats1_up, bt_hfstats1_down, bt_hfstats2_up, bt_hfstats2_down  # 22-25
                ]
        
    @staticmethod
    def ttbb_worker(roo):
        with uproot.open(roo) as f:
            t = f['Events']
            gw    = t.array('genWeight'  )
            ttbb  = t.array('genTtbarId')
            ttbb  = ttbb % 100
            scale = t.array('LHEScaleWeight').pad(9).fillna(1)
            ps    = t.array('PSWeight'      ).pad(4).fillna(1)
            # need pdf as well
            pdf_up   = t.array('pdfWeight_Up') * np.sign(gw)
            pdf_down = t.array('pdfWeight_Down') * np.sign(gw)
            tot_count= sum(gw>=0.0) - sum(gw<0.0)
            #
            bt_nom = sum(t.array('BTagWeight') * np.sign(gw))
            bt_jes_up, bt_jes_down = sum(t.array('BTagWeight_jes_up') * np.sign(gw)), sum(t.array('BTagWeight_jes_down') * np.sign(gw))
            bt_lf_up, bt_lf_down =   sum(t.array('BTagWeight_lf_up') * np.sign(gw)),  sum(t.array('BTagWeight_lf_down') * np.sign(gw))
            bt_hf_up, bt_hf_down =   sum(t.array('BTagWeight_hf_up') * np.sign(gw)),  sum(t.array('BTagWeight_hf_down') * np.sign(gw))
            bt_lfstats1_up, bt_lfstats1_down = sum(t.array('BTagWeight_lfstats1_up') * np.sign(gw)), sum(t.array('BTagWeight_lfstats1_down') * np.sign(gw))
            bt_lfstats2_up, bt_lfstats2_down = sum(t.array('BTagWeight_lfstats2_up') * np.sign(gw)), sum(t.array('BTagWeight_lfstats2_down') * np.sign(gw))
            bt_hfstats1_up, bt_hfstats1_down = sum(t.array('BTagWeight_hfstats1_up') * np.sign(gw)), sum(t.array('BTagWeight_hfstats1_down') * np.sign(gw))
            bt_hfstats2_up, bt_hfstats2_down = sum(t.array('BTagWeight_hfstats2_up') * np.sign(gw)), sum(t.array('BTagWeight_hfstats2_down') * np.sign(gw))
            #
            ttbb_count = sum((gw>=0.0) & (ttbb>=51)) - sum((gw<0.0) & (ttbb>=51))
            #
            scale = np.sign(gw) * scale
            sc_tot = [ sum( scale[:,i] ) for i in range(9)]
            mur_up_tot, mur_down_tot   = sc_tot[7], sc_tot[1]
            muf_up_tot, muf_down_tot   = sc_tot[5], sc_tot[3]
            murf_up_tot, murf_down_tot = sc_tot[8], sc_tot[0]
            sc_ttbb = [ sum( scale[(ttbb>=51)][:,i] ) for i in range(9)]
            mur_up_ttbb, mur_down_ttbb   = sc_ttbb[7], sc_ttbb[1]
            muf_up_ttbb, muf_down_ttbb   = sc_ttbb[5], sc_ttbb[3]
            murf_up_ttbb, murf_down_ttbb = sc_ttbb[8], sc_ttbb[0]
            #
            ps = np.sign(gw) * ps
            ps_tot = [ sum( ps[:,i] ) for i in range(4) ]
            isr_up_tot, isr_down_tot = ps_tot[2], ps_tot[0]
            fsr_up_tot, fsr_down_tot = ps_tot[3], ps_tot[1]
            ps_ttbb = [ sum( ps[(ttbb>=51)][:,i] ) for i in range(4) ]
            isr_up_ttbb, isr_down_ttbb = ps_ttbb[2], ps_ttbb[0]
            fsr_up_ttbb, fsr_down_ttbb = ps_ttbb[3], ps_ttbb[1]
            #
            pdf_up_tot, pdf_down_tot = sum(pdf_up), sum(pdf_down)
            pdf_up_ttbb, pdf_down_ttbb = sum(pdf_up[(ttbb>=51)]), sum(pdf_down[(ttbb>=51)])
            #
            return [tot_count,  # 0
                    mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot, # 1-6
                    isr_up_tot, isr_down_tot, fsr_up_tot, fsr_down_tot, # 7-10
                    pdf_up_tot, pdf_down_tot, # 11-12
                    ttbb_count,  # 13
                    mur_up_ttbb, mur_down_ttbb, muf_up_ttbb, muf_down_ttbb, murf_up_ttbb, murf_down_ttbb, # 14-19
                    isr_up_ttbb, isr_down_ttbb, fsr_up_ttbb, fsr_down_ttbb,# 20-23
                    pdf_up_ttbb, pdf_down_ttbb, # 24-25
                    bt_nom, # 26
                    bt_jes_up, bt_jes_down, bt_lf_up, bt_lf_down, bt_hf_up, bt_hf_down, # 27-32
                    bt_lfstats1_up, bt_lfstats1_down, bt_lfstats2_up, bt_lfstats2_down,  # 33-36
                    bt_hfstats1_up, bt_hfstats1_down, bt_hfstats2_up, bt_hfstats2_down # 37-40
                ]

if __name__ == '__main__':
    test = PreSkim('WWW','2016') # sample, year, isttbar=False, isttbb=False
    print(test.get_metadata())
    
