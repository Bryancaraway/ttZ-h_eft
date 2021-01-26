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
from modules.PostSkim import PostSkim
from lib.fun_library import clop_pear_ci

eps = 0.0001

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
    sf_period_dict = {
        'Electron':{'2016':['2016_AtoF','2016_GtoH'],
                    '2017':['2017'],
                    '2018':['2018_before','2018_after']},
        'Muon':{'2016':['2016'],
                '2017':['2017'],
                '2018':['2018']}
    }
        
                    

    def __init__(self,channel, year):
        self.channel = channel # Electron or Muon
        self.year = year
        # load data and mc that passes reference trigger
        self.data_df = self.getData('Data_Single'+self.of_lep[self.channel])
        di_mc_df, di_gen_df   = self.getData('TTTo2L2Nu')
        semi_mc_df, semi_gen_df   = self.getData('TTToSemiLeptonic')
        self.mc_df  = pd.concat([di_mc_df,semi_mc_df], axis='rows', ignore_index=True)
        self.gen_df = AnaDict({k: PostSkim.try_concatenate([di_gen_df[k],semi_gen_df[k]]) for k in di_gen_df})
        # need to load and add sf, topptrwgt to mc_df
        for ref_trig_period in self.sf_period_dict[self.channel][self.year]:
            self.add_refTrigSF(ref_trig_period)
        self.add_lepSF()
        self.add_toppt_reweight()
        self.get_total_weight()
        self.calc_trig_eff_sf()
        #self.plot_test('nBottoms',[-.5,0.5,1.5,2.5,3.5,4.5])
        #self.plot_test('Muon_pt',[30,  40,  50,  60, 120, 200, 500])
        # self.plot_test('Muon_eta',[0,  0.9, 1.2, 2.1, 2.4])
        #pt_bins  = [30,  40,  50,  60, 120, 200, 500]
        #eta_bins = [0,  0.9, 1.2, 2.1, 2.4]


    def calc_trig_eff_sf(self):
        bins_dict = {
            'Electron': {'pt': [30,  40,  55, 120, 200, 500], 'eta': [-2.5,-1.5,-0.8,0.,0.8,1.5,2.5]},
            'Muon'    : {'pt':[30,  40,  50,  60, 120, 200, 500], 'eta': [-2.4,-2.1,-1.2,-0.9,0.0,0.9,1.2,2.1,2.4]}
        }
        bins = bins_dict[self.channel]
        num_mc, edges = np.histogram(self.mc_df[self.channel+'_pt'][cfg.hlt_path[self.channel.lower()][self.year](self.mc_df)].clip(bins['pt'][0], bins['pt'][-1]), 
                                     bins=bins['pt'], weights=self.mc_df['tot_weight'][cfg.hlt_path[self.channel.lower()][self.year](self.mc_df)])
        den_mc,   _   = np.histogram(self.mc_df[self.channel+'_pt'].clip(bins['pt'][0], bins['pt'][-1]), bins=bins['pt'], weights=self.mc_df['tot_weight'])
        #
        num_data, _   = np.histogram(self.data_df[self.channel+'_pt'][cfg.hlt_path[self.channel.lower()][self.year](self.data_df)].clip(bins['pt'][0], bins['pt'][-1]), 
                                     bins=bins['pt'])
        den_data, _   = np.histogram(self.data_df[self.channel+'_pt'].clip(bins['pt'][0], bins['pt'][-1]), bins=bins['pt'])
        print("MC ratio per pt bin",num_mc/den_mc)
        print("Data ratio per pt bin",num_data/den_data)



    def plot_test(self,k_name,bins):
        n_mc, edges = np.histogram(self.mc_df[  k_name].clip(bins[0], bins[-1]) , bins=bins, weights=self.mc_df['tot_weight'])
        n_mc2,  _   = np.histogram(self.mc_df[  k_name].clip(bins[0], bins[-1]) , bins=bins, weights=self.mc_df['tot_weight']**2)
        n_data, _   = np.histogram(self.data_df[k_name].clip(bins[0], bins[-1]) , bins=bins, )
        yerr = (n_data/n_mc)*np.sqrt((np.sqrt(n_data)/n_data)**2+ (np.sqrt(n_mc2)/n_mc)**2)
        plt.errorbar(x=(edges[1:]+edges[:-1])/2, y=n_data/n_mc, xerr=(edges[1:]-edges[:-1])/2, yerr=yerr, fmt='.', color='k')
        plt.errorbar(x=(edges[1:]+edges[:-1])/2, y=n_mc/n_mc, xerr=(edges[1:]-edges[:-1])/2, yerr=np.sqrt(n_mc2)/n_mc, fmt='.', color='red')
        plt.axhline(1, color='red', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        plt.xlim(bins[0],bins[-1])
        plt.ylabel('Data/MC')
        plt.ylim(0.5,1.5)
        plt.xlabel(k_name)
        plt.title('Data/MC using only ref trigger (stat err. only)')
        plt.show()
        

    def add_refTrigSF(self,period):
        file_sf_dict = {
            'Electron':{'2016':(lambda : json.load(open(".json"))),
                        '2017':(lambda : json.load(open(".json"))),
                        '2018':(lambda : json.load(open(".json"))),
                    },
            'Muon':    {
                '2016_AtoF'  : (lambda : uproot.open(
                    f"{cfg.dataDir}/lep_effsf_files/EfficienciesAndSF_RunBtoF.root")['IsoMu24_OR_IsoTkMu24_PtEtaBins']['pt_abseta_ratio']),
                '2016_GtoH'  : (lambda : uproot.open(
                    f"{cfg.dataDir}/lep_effsf_files/EfficienciesAndSF_Period4.root")['IsoMu24_OR_IsoTkMu24_PtEtaBins']['pt_abseta_ratio']),
                '2017'       : (lambda : uproot.open(
                    f"{cfg.dataDir}/lep_effsf_files/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root")['IsoMu27_PtEtaBins']['pt_abseta_ratio']),
                '2018_before': (lambda : uproot.open(
                    f"{cfg.dataDir}/lep_effsf_files/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root")['IsoMu24_PtEtaBins']['pt_abseta_ratio']),
                '2018_after' : (lambda : uproot.open(
                    f"{cfg.dataDir}/lep_effsf_files/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root")['IsoMu24_PtEtaBins']['pt_abseta_ratio']),
            }
        }
        self.ref_trigger = file_sf_dict[self.of_lep[self.channel]][period]()
        if self.channel == 'Electron':
            self.handleMu_refTrigSF(period)
        else:
            self.handleEle_refTrigSF()
            
    def handleMu_refTrigSF(self,period):
        rel_lumi_period = {
            '2016_AtoF':1-0.4518,
            '2016_GtoH':.4518,
            '2017':1,
            '2018_before':76/477,
            '2018_after' :401/477,
        }
        lep_pt, lep_eta = self.mc_df['Muon_pt'], abs(self.mc_df['Muon_eta'])
        hist = self.ref_trigger
        pt_bins, eta_bins = hist.edges
        v    = np.array(hist.values).T
        pt_digi  = pd.cut(lep_pt.clip(pt_bins[0]+eps,pt_bins[-1]-eps), bins=pt_bins, right=True, include_lowest=True, labels=range(len(pt_bins)-1))
        eta_digi = pd.cut(np.clip(lep_eta,eta_bins[0]+eps,eta_bins[-1]-eps), bins=eta_bins, right=True, include_lowest=True, labels=range(len(eta_bins)-1))
        if 'refTrig_sf' in self.mc_df:
            self.mc_df['refTrig_sf'] += (np.array([ v[x][y] for x,y in zip(eta_digi.values,pt_digi.values)]) * rel_lumi_period[period] )
        else:
            self.mc_df['refTrig_sf'] = (np.array([ v[x][y] for x,y in zip(eta_digi.values,pt_digi.values)]) * rel_lumi_period[period] )

    def handleEle_refTrigSF(self):
        pass
        

    def getData(self,sample):
        input_file = cfg.postSkim_dir+f"{self.year}/Trig_{sample_cfg[sample]['out_name']}/{sample}.pkl"
        isData   = 'Data' in sample
        gD_out = AnaDict.read_pickle(f'{input_file}')
        pass_ref = self.pass_refTrig(gD_out['events'])
        pass_nb  = (lambda df: (df['nBottoms'] == 1) | (df['nBottoms'] == 2))(gD_out['events'])
        pass_event = pass_ref & pass_nb
        if isData:
            return gD_out['events'][pass_event]
        else:
            events = gD_out['events'][pass_event]
            gen    = AnaDict(gD_out['gen'])[pass_event.values]
            return (events, gen)
        #if isData:
        #    self.data_df = pd.concat([self.data_df,gD_out['events']], axis='rows', ignore_index=True)
    def pass_refTrig(self,df):
        ref_trig_dict = {
            'Electron':{'2016': (lambda x : (x['HLT_IsoMu24'] | 
                                             x['HLT_IsoTkMu24'])),
                        '2017': (lambda x : (x['HLT_IsoMu27'])),
                        '2018': (lambda x : (x['HLT_IsoMu24'])),
                    },
            'Muon'    : cfg.hlt_path['electron'],
        }
        pass_reftrig = ref_trig_dict[self.channel][self.year](df)
        return pass_reftrig
        
    def add_lepSF(self):
        for lep, lep_type in zip(['Muon','Electron'],['muon','electron']):
            lep_pt  = self.mc_df[f'{lep}_pt']
            lep_eta = self.mc_df[f'{lep}_eta'] 
            sf = AnaDict.read_pickle(f'{cfg.dataDir}lep_sf_files/{lep_type}_sf_{self.year}.pkl')
            total_lep_sf = [np.ones(len(lep_pt)),np.zeros(len(lep_pt))]
            for k in sf:
                print(k)
                # dict structure is sf type : eta_bins: pt_bins: values, up, down
                if ((lep_type == 'muon' and self.year != '2016') or k == 'SF'): lep_eta = abs(lep_eta)
                eta_keys = list(sf[k].keys())
                pt_keys  = list(sf[k][eta_keys[-1]].keys())
                eta_bins=[float(bin_str.split(',')[0]) for bin_str in eta_keys] + [float(eta_keys[-1].split(',')[1])]
                pt_bins =[float(bin_str.split(',')[0]) for bin_str in pt_keys] +  [float( pt_keys[-1].split(',')[1])]
                eta_digi = pd.cut(np.clip(lep_eta,eta_bins[0]+eps,eta_bins[-1]-eps), bins=eta_bins, right=True, include_lowest=True,
                                  labels=sorted(eta_keys, key= lambda z : float(z.split(',')[0])))
                pt_digi  = pd.cut(lep_pt.clip(pt_bins[0]+eps,pt_bins[-1]-eps), bins=pt_bins, right=True, include_lowest=True, 
                                  labels=sorted(pt_keys, key= lambda z : float(z.split(',')[0])))
                lep_sf = np.array([(lambda x,y: np.array([sf[k][x][y]['values'], sf[k][x][y]['up']]))(x,y) for x,y in zip(eta_digi.values,pt_digi.values)]) 
                total_lep_sf[0] *= lep_sf[:,0]
                lep_sf_err = abs(lep_sf[:,1] - lep_sf[:,0]) 
                total_lep_sf[1] = total_lep_sf[0]*np.sqrt(np.power(lep_sf_err/lep_sf[:,0],2)+ np.power(total_lep_sf[1]/total_lep_sf[0],2)) 
            self.mc_df[f'{lep_type}_sf']      = total_lep_sf[0]
            self.mc_df[f'{lep_type}_sf_up']   = total_lep_sf[1]+total_lep_sf[0]
            self.mc_df[f'{lep_type}_sf_down'] = total_lep_sf[0]-total_lep_sf[1]

    def add_toppt_reweight(self):
        tt_pt = self.gen_df['GenPart_pt'][(abs(self.gen_df['GenPart_pdgId']) == 6)]
        # Using the newer theo (NNLO QCD + NLO EW) corrections which is better for BSM analysis aspects
        sf = (lambda x: 0.103*np.exp(-0.0118*np.clip(x,0,np.inf)) - 0.000134*np.clip(x,0,np.inf) + 0.973) 
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Case_3_3_The_Effective_Field_The
        #https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf'
        # the theo toppt event re-weighting unc. is based on [1, w**2] where w is the event reweighting 
        toppt_rwgt = np.sqrt(sf(tt_pt[:,0]) * sf(tt_pt[:,1])) 
        toppt_rwgt_up = np.where(toppt_rwgt > 1.0, toppt_rwgt**2,  1.0)
        toppt_rwgt_dn = np.where(toppt_rwgt < 1.0, toppt_rwgt**2,  1.0)
        self.mc_df['topptWeight']      = toppt_rwgt
        self.mc_df['topptWeight_Up']   = toppt_rwgt_up
        self.mc_df['topptWeight_Down'] = toppt_rwgt_dn

    def get_total_weight(self):
        self.mc_df['tot_weight'] = (self.mc_df['weight']* np.sign(self.mc_df['genWeight']) 
                                    #* (np.where(self.mc_df['weight']>300,0,1))
                                    * self.mc_df['topptWeight']
                                    * (self.mc_df['HEM_weight']  if self.year == '2018' else 1.0 )
                                    * self.mc_df['muon_sf']
                                    * self.mc_df['electron_sf']
                                    * self.mc_df['refTrig_sf']
                                    * self.mc_df['BTagWeight'] 
                                    * self.mc_df['puWeight']  
                                    * (self.mc_df['PrefireWeight'] if self.year != '2018' else 1.0))
        

if __name__ == '__main__':
    Calc_LepEffSF('Electron','2016')
    Calc_LepEffSF('Electron','2017')
    Calc_LepEffSF('Electron','2018')
