########################
### Skim metadata    ###
### for analysis     ###
########################
### written by:      ###
### Bryan Caraway    ###
########################
##
#
from uproot_methods import TLorentzVectorArray
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

    def __init__(self, sample, year, isData, tree, jec_sys=None, run_norm=True):
        self.sample = sample
        self.year   = year
        self.isData = isData
        self.tree   = tree
        #
        self.isttbar = 'TTTo' in self.sample or 'TTJets' in self.sample
        self.isttbb  = 'TTbb' in self.sample
        self.issig   = 'TTZ' in self.sample or 'ttH' in self.sample or 'TTH' in self.sample
        #self.metadata = {'sample': self.sample, 'year': self.year, 
        #                 'xs': sample_cfg[self.sample]['xs'], 'kf': sample_cfg[self.sample]['kf']}
        self.metadata = {}
        #
        if run_norm:
            if not self.isData:
                if jec_sys is None:
                    if self.isttbar or self.isttbb:
                        self.ttbar_ttbb_normalization()
                    elif self.issig:
                        self.event_signal_normalization()
                    else:
                        self.event_var_normalization()
                else:
                    self.jec_var_normalization()

    ## main functions    
    
    def get_metadata(self):
        return self.metadata

    def jec_var_normalization(self):
        results = np.array(self.jec_worker(self.tree))
        self.metadata.update({'tot_events':results[0]})

    def event_var_normalization(self): # 1,2,4
        results = np.array(self.event_worker(self.tree))
        self.metadata.update({
            'tot_events' :results[0], 
            'mur_up_tot' :results[1], 'mur_down_tot' :results[2],
            'muf_up_tot' :results[3], 'muf_down_tot' :results[4],
            'murf_up_tot':results[5], 'murf_down_tot':results[6],
            'pdf_up_tot' :results[7], 'pdf_down_tot' :results[8],
        })
    def event_signal_normalization(self): # 1,2,4
        results = np.array(self.signal_worker(self.tree, self.sample))
        genpt_bins = [0,200,300,450,np.inf]
        keys = ['events_tot','events_rap_tot',
                'mur_up_tot','mur_down_tot','muf_up_tot','muf_down_tot','murf_up_tot','murf_down_tot',
                'isr_up_tot','isr_down_tot','fsr_up_tot','fsr_down_tot',
                'pdf_up_tot','pdf_down_tot']
        _upd_dict = {}
        for i,ptbin in enumerate(genpt_bins):
            s_format = (lambda x: x.replace('_tot',f'_{i-1}')) if i > 0 else (lambda x:x)
            for j,key in enumerate(keys):
                _upd_dict[s_format(key)] = results[i*len(keys)+j]
        #
        self.metadata.update({'tot_events':results[0]})
        self.metadata.update( _upd_dict )

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
    def jec_worker(t):
        gw = t.array('genWeight')
        tot_count= sum(gw>=0.0) - sum(gw<0.0)
        return [tot_count]

    @staticmethod
    def event_worker(t):
        scps_helper = (lambda arr: np.clip(np.nanpercentile(arr,.15), np.nanpercentile(arr,99.85), arr))
        gw = t.array('genWeight')#, executor=executor, blocking=True)
        try:
            scale = t.array('LHEScaleWeight').pad(9).fillna(1) * np.sign(gw)
        except:
            scale = np.ones(shape=(len(gw),9))
        pdf_up   = t.array('pdfWeight_Up') * np.sign(gw)
        pdf_down = t.array('pdfWeight_Down') * np.sign(gw)
        #
        sc_tot = [ sum( scps_helper(scale[:,i] )) for i in range(9)]
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
    def signal_worker(t, sample):
        scps_helper = (lambda arr: np.clip(np.nanpercentile(arr,.15), np.nanpercentile(arr,99.85), arr))
        genpt_bins = [0,200,300,450,np.inf] # 500 is overflow
        gen_id, gen_st, gen_mom = [t.array(k) for k in ['GenPart_pdgId','GenPart_status','GenPart_genPartIdxMother']]
        zhcut = ( ((gen_id == 23)|(gen_id == 25)) & (gen_st == 62))
        id_cut = (zhcut.sum() > 0).flatten()
        if sample == 'TTZToQQ':
            isqq_fromZ = (((abs(gen_id) <  5) & (gen_id[gen_mom] == 23) & (gen_st[gen_mom] == 62)).sum() == 2).flatten()
            id_cut = ((id_cut==True) & (isqq_fromZ==True))
        #print(len(zhcut), len(id_cut), len(zhcut[id_cut]))
        zhcut = zhcut[id_cut]
        tarray = (lambda k : t.array(k)[id_cut])
        del gen_id, gen_st
        gw = tarray('genWeight')#, executor=executor, blocking=True)
        zh_pt, zh_eta, zh_phi, zh_mass = [tarray(k)[zhcut].flatten() for k in ['GenPart_pt','GenPart_eta','GenPart_phi','GenPart_mass']]
        getTLVm = TLorentzVectorArray.from_ptetaphim
        zh_tlv   = getTLVm(zh_pt,zh_eta,zh_phi,zh_mass)
        stxs_cut = abs(zh_tlv.rapidity) <= 2.5
        del zh_tlv, zh_eta, zh_phi, zh_mass
        #
        pdf_up   = tarray('pdfWeight_Up') * np.sign(gw)
        pdf_down = tarray('pdfWeight_Down') * np.sign(gw)
        try:
            scale = tarray('LHEScaleWeight').pad(9).fillna(1)
        except:
            scale = np.ones(shape=(len(gw),9))
        try:
            ps = tarray('PSWeight').pad(4).fillna(1)
        except:
            ps = np.ones(shape=(len(gw),4))
        #
        #ps = ps * np.sign(gw)
        #scale = scale * np.sign(gw)
        ps    = np.array([scps_helper(ps[:,i]) * np.sign(gw) for i in range(4)]).T
        scale = np.array([scps_helper(scale[:,i]) * np.sign(gw) for i in range(9)]).T
        #print(ps.shape)
        out_list = []
        for i,ptbin in enumerate(genpt_bins):
            if ptbin == 0: # do inclusive counts first 
                pt_cut = (zh_pt >= 0)
            else:
                pt_cut = ((zh_pt >= genpt_bins[i-1]) & (zh_pt < ptbin))
            #
            #print(len(gw),len(pt_cut), len(zh_pt), len(zhcut), len(id_cut))
            #65459 1705123 1705123 65459 71774
            genpt_tot = sum(gw[pt_cut==True]>=0.0) - sum(gw[pt_cut==True]<0.0) 
            genpt_rap_tot = sum(gw[(pt_cut==True) & (stxs_cut==True)]>=0.0) - sum(gw[(pt_cut==True) & (stxs_cut==True)]<0.0) 
            ps_tot = [ sum( ps[:,i][(pt_cut==True) & (stxs_cut==True)] ) for i in range(4) ]
            #ps_tot = [ sum( scps_helper(ps[:,i][(pt_cut==True) & (stxs_cut==True)]) ) for i in range(4) ]
            isr_up_tot, isr_down_tot = ps_tot[2], ps_tot[0]
            fsr_up_tot, fsr_down_tot = ps_tot[3], ps_tot[1]
            #[0] is ISR=0.5 FSR=1; [1] is ISR=1 FSR=0.5; [2] is ISR=2 FSR=1; [3] is ISR=1 FSR=2
            #
            #sc_tot = [ sum( scps_helper(scale[:,i][(pt_cut==True) & (stxs_cut==True)]) ) for i in range(9)]
            sc_tot = [ sum( scale[:,i][(pt_cut==True) & (stxs_cut==True)] ) for i in range(9)]
            mur_up_tot, mur_down_tot   = sc_tot[7], sc_tot[1]
            muf_up_tot, muf_down_tot   = sc_tot[5], sc_tot[3]
            murf_up_tot, murf_down_tot = sc_tot[8], sc_tot[0]
            #
            pdf_up_tot, pdf_down_tot = map((lambda x : sum(x[(pt_cut==True) & (stxs_cut==True)])), [pdf_up,pdf_down])
            #
            out_list += [genpt_tot, genpt_rap_tot, 
                         mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot,
                         isr_up_tot, isr_down_tot, fsr_up_tot, fsr_down_tot,
                         pdf_up_tot, pdf_down_tot] 
        #
        return out_list
        
    @staticmethod
    def ttbb_worker(t):
        scps_helper = (lambda arr: np.clip(np.nanpercentile(arr,.15), np.nanpercentile(arr,99.85), arr))
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
        sc_tot = [ sum( scps_helper(scale[:,i] )) for i in range(9)]
        mur_up_tot, mur_down_tot   = sc_tot[7], sc_tot[1]
        muf_up_tot, muf_down_tot   = sc_tot[5], sc_tot[3]
        murf_up_tot, murf_down_tot = sc_tot[8], sc_tot[0]
        sc_ttbb = [ sum( scps_helper(scale[(ttbb>=51)][:,i] )) for i in range(9)]
        mur_up_ttbb, mur_down_ttbb   = sc_ttbb[7], sc_ttbb[1]
        muf_up_ttbb, muf_down_ttbb   = sc_ttbb[5], sc_ttbb[3]
        murf_up_ttbb, murf_down_ttbb = sc_ttbb[8], sc_ttbb[0]
        #
        ps = np.sign(gw) * ps
        ps_tot = [ sum( scps_helper(ps[:,i] )) for i in range(4) ]
        isr_up_tot, isr_down_tot = ps_tot[2], ps_tot[0]
        fsr_up_tot, fsr_down_tot = ps_tot[3], ps_tot[1]
        ps_ttbb = [ sum( scps_helper(ps[(ttbb>=51)][:,i] )) for i in range(4) ]
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
            active_mask = (lambda df : df[p_mask]) if p_mask is not None else (lambda df : df)
            opt         = f'{opt}_' if opt is not None else ''
            n_ak4jets = np.clip(0, 9, active_mask(events['n_ak4jets']))
            ht        = np.clip(0,1500, active_mask(jets['Jet_pt']).sum())
            b_score   = np.clip(0,1, active_mask(jets['Jet_btagDeepB'])[:,0])
            gw_np     = active_mask(np.sign(events['genWeight']))
            #print(n_ak4jets.shape, ht.shape, b_score.shape, gw_np.shape)
            n,_ = np.histogram(n_ak4jets, bins=np.arange(-0.5,10.5,1), weights=gw_np)
            ht_vs_n,*_ = np.histogram2d(x=n_ak4jets, y=ht, bins=[np.arange(-0.5,10.5,1),np.arange(0,1550,50)], weights=gw_np)
            b_vs_n, *_ = np.histogram2d(x=n_ak4jets, y=b_score, bins=[np.arange(-0.5,10.5,1),np.arange(0,1.05,0.05)], weights=gw_np)
            self.metadata[f'{opt}nj_yield'] = n
            # for closure
            self.metadata[f'{opt}htvsn_2d'] = ht_vs_n 
            self.metadata[f'{opt}bvsn_2d'] = b_vs_n
            ht_vs_n_bw,*_ = np.histogram2d(x=n_ak4jets, y=ht, bins=[np.arange(-0.5,10.5,1),np.arange(0,1550,50)], 
                                          weights=active_mask(jets['Jet_btagSF_deepcsv_shape'].prod())*gw_np)
            b_vs_n_bw, *_ = np.histogram2d(x=n_ak4jets, y=b_score, bins=[np.arange(-0.5,10.5,1),np.arange(0,1.05,0.05)], 
                                          weights=active_mask(jets['Jet_btagSF_deepcsv_shape'].prod())*gw_np)
            self.metadata[f'{opt}htvsn_2dsf'] = ht_vs_n_bw
            self.metadata[f'{opt}bvsn_2dsf'] = b_vs_n_bw
            # 
            jetshape_sf_names = re.findall(r'Jet_btagSF_deepcsv_shape\w*' ,' '.join(jets.keys()))#Jet_btagSF_deepcsv_shape_up_cferr2
            #btag_w_names = []
            for sf_n in jetshape_sf_names:
                btag_w_name = re.search(r'shape\w*', sf_n).group().replace('shape','BTagWeight')
                #btag_w_names.append(btag_w_name)
                events[btag_w_name] = jets[sf_n].prod() # add weight to events
                #for btag_w_name in btag_w_names: # get yields and add to metadata
                nbw,_ = np.histogram(n_ak4jets, bins=np.arange(-0.5,10.5,1), 
                                   weights=(active_mask(events[btag_w_name])*gw_np))
                self.metadata[f'{opt}btw_yield'+btag_w_name.replace('BTagWeight','')] = nbw
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
