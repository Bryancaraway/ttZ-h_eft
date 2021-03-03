########################
### Skim metadata    ###
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
from config.sample_cff import sample_cfg
#
import numpy as np
np.random.seed(0)
import pandas as pd


class SkimMeta :
    '''
    only run full code on MC
    get 4 main items from sample
    1) event normalization  (all) 
    2) varied normalization (all)
    3) varied normalization w/resp tt+bb (only ttbar and ttbb)
    4) btagSFweight yield per np, and per jet multi
    '''
    #outDir = '' # add later */*postSkim
    #nanoAODv7_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/'
    #sample  = None

    def __init__(self, sample, year, isData, tree):
        self.sample = sample
        self.year   = year
        self.isData = isData
        self.tree   = tree
        #
        self.isttbar = 'TTTo' in self.sample or 'TTJets' in self.sample
        self.isttbb  = 'TTbb' in self.sample
        #self.metadata = {'sample': self.sample, 'year': self.year, 
        #                 'xs': sample_cfg[self.sample]['xs'], 'kf': sample_cfg[self.sample]['kf']}
        self.metadata = {}
        #
        if not self.isData:
            if self.isttbar or self.isttbb:
                self.ttbar_ttbb_normalization()
            else:
                self.event_var_normalization()

    ## main functions    
    
    def get_metadata(self):
        return self.metadata

    def event_var_normalization(self): # 1,2,4
        results = np.array(self.event_worker(self.tree))
        self.metadata.update({
            'tot_events' :results[0], 
            'mur_up_tot' :results[1], 'mur_down_tot' :results[2],
            'muf_up_tot' :results[3], 'muf_down_tot' :results[4],
            'murf_up_tot':results[5], 'murf_down_tot':results[6],
            'pdf_up_tot' :results[7], 'pdf_down_tot' :results[8],
        })

    def ttbar_ttbb_normalization(self): # 3,4
        results = np.array(self.ttbb_worker(self.tree))
        self.metadata.update({
            'tot_events' :results[0], 
            'mur_up_tot' :results[1], 'mur_down_tot'  :results[2],
            'muf_up_tot' :results[3], 'muf_down_tot'  :results[4],
            'murf_up_tot':results[5], 'murf_down_tot' :results[6],
            'isr_up_tot' :results[7], 'isr_down_tot'  :results[8],
            'fsr_up_tot' :results[9], 'fsr_down_tot'  :results[10],
            'pdf_up_tot' :results[11],'pdf_down_tot'  :results[12],
            #
            'ttbb_events' : results[13], 
            'mur_up_ttbb' :results[14],'mur_down_ttbb' :results[15],
            'muf_up_ttbb' :results[16],'muf_down_ttbb' :results[17],
            'murf_up_ttbb':results[18],'murf_down_ttbb':results[19],
            'isr_up_ttbb' :results[20],'isr_down_ttbb' :results[21],
            'fsr_up_ttbb' :results[22],'fsr_down_ttbb' :results[23],
            'pdf_up_ttbb' :results[24],'pdf_down_ttbb' :results[25],
        })
    

    ## worker functions 
    @staticmethod
    def event_worker(t):
        gw = t.array('genWeight')#, executor=executor, blocking=True)
        try:
            scale = t.array('LHEScaleWeight').pad(9).fillna(1) * np.sign(gw)
        except:
            scale = np.ones(shape=(len(gw),9))
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
        tot_count= sum(gw>=0.0) - sum(gw<0.0)
        #
        return [tot_count,        # 0
                mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot, # 1-6
                pdf_up_tot, pdf_down_tot, # 7-8
        ]
        
    @staticmethod
    def ttbb_worker(t):
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
        ]

    def add_btagweightsf_counts(self, jets, events, geninfo):
        # need to only count certain processes/sub-processes for later
        # main things to worry about tt+bb, tt+lf, ttZbb in TTZToQQ
        ## -- ## 
        def btagweightsf_helper(p_mask=None, opt=None):
            active_mask = (lambda df : df[p_mask]) is p_mask is not None else (lambda df : df)
            opt         = f'{opt}_' if opt is not None else ''
            n_ak4jets = np.clip(0, 12, active_mask(events['n_ak4jets']))
            gw_np     = active_mask(np.sign(events['genWeight']))
            n,_ = np.histogram(n_ak4jets, bins=np.arange(-0.5,13.5,1), weights=gw_np)
            self.metadata[f'{opt}nj_yield'] = n
            jetshape_sf_names = re.findall(r'Jet_btagSF_deepcsv_shape\w*' ,' '.join(jets.keys()))#Jet_btagSF_deepcsv_shape_up_cferr2
            #btag_w_names = []
            for sf_n in jetshape_sf_names:
                btag_w_name = re.search(r'shape\w*', sf_n).group().replace('shape','BTagWeight')
                #btag_w_names.append(btag_w_name)
                events[btag_w_name] = jets[sf_n].prod() # add weight to events
                #for btag_w_name in btag_w_names: # get yields and add to metadata
                nbw,_ = np.histogram(n_ak4jets, bins=np.arange(-0.5,13.5,1), 
                                   weights=(active_mask(events[btag_w_name])*gw_np))
                self.metadata['{opt}btw_yield'+btag_w_name.replace('BTagWeight','')] = nbw
                # no longer need shape sf
                del jets[sf_n]
        ## -- ##
        ##### may have to add special case for tt+cc from ttjets (5FS) sample #####
        process_mask = None
        if   'TTbb'    in self.sample:
            process_mask = (geninfo['genTtbarId'] % 100 >= 51)
            btagweightsf_helper(process_mask)
        elif 'TTTo'    in self.sample:
            process_mask = (geninfo['genTtbarId'] % 100 < 51)
            # process_mask = [(geninfo['genTtbarId'] % 100 < 41),((geninfo['genTtbarId'] % 100 < 50) & (geninfo['genTtbarId'] % 100 >= 41))]
            # for p_, mask_ in zip(['lf','cc'],process_mask)
            #     btagweightsf_helper(mask_, p_)
            btagweightsf_helper(process_mask)
        elif 'TTZToQQ' in self.sample:
            isbb_fromZ     = ((abs(geninfo['GenPart_pdgId']) == 5) & (geninfo['GenPart_pdgId'][geninfo['GenPart_genPartIdxMother']] == 23))
            isZbb  = ((geninfo['GenPart_pdgId'] == 23) & (isbb_fromZ.sum() == 2))
            process_mask = isZbb.sum() == 0
            btagweightsf_helper(process_mask)
        else:
            btagweightsf_helper()
        

if __name__ == '__main__':
    test = PreSkim('WWW','2016') # sample, year, isttbar=False, isttbb=False
    print(test.get_metadata())
