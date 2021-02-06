import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import uproot
import numpy as np
import pandas as pd
import json
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
from modules.AnaDict import AnaDict
from lib.fun_library import save_pdf, clop_pear_ci
import config.ana_cff as cfg

# use plotter class to load all necessary data
class DeepAk8_fakeSF:
    bbvl_sf_file = cfg.dataDir+'/deepak8sf/deepak8_bbvL_sf.json'
    
    def __init__(self,year):
        self.year = year
        self.sf_dict = {'sig_sf':self.get_sf_from_json()}
        Plotter.load_data(year, samples=cfg.Sig_MC+cfg.Bkg_MC, byprocess=True)
        self.mc_dict   = Plotter.data_dict
        # seperate between mc 1 and mc 2 
        self.initial_datamc_cuts = (lambda df_:( ( (df_['isEleE']==True) | (df_['isMuonE']==True) ) &
                                                 (df_['Zh_bbvLscore'] >= 0.0) & (df_['n_b_outZh']>=2)))
        self.mc1_df, self.mc2_df = self.seperate_mc_dict()
        self.data_df = self.get_data_df()
        # get f_sf
        self.get_f_sf()
        self.store_bbvl_sf()
    

    def get_f_sf(self):
        #bbvl_cut = (lambda _df: 
        mc1_df_copy = self.mc1_df.copy()
        self.mc1_df.loc[:,'tot_weight'] = self.mc1_df['tot_weight']*self.mc1_df['dak8md_bbvl_sf']
        self.data_df['tot_weight'] = 1
        #
        def get_eff(_df):
            _num, *edges = np.histogram2d(x= _df['Zh_bbvLscore'], y= _df['Zh_pt'].clip(self.pt_bins[0],self.pt_bins[-1]), bins=[[0.]+self.score_bins,self.pt_bins], weights= _df['tot_weight'])
            _den  = np.nansum(_df['tot_weight'])
            _eff_dwn, _eff_up = clop_pear_ci(_num,_den)
            return _num, _den, _eff_up, _eff_dwn, edges
        #
        num_mc1sf, den_mc1sf, *_ = get_eff(self.mc1_df)
        num_mc1, den_mc1, *_ = get_eff(mc1_df_copy)
        real_sf1 = (den_mc1 - np.sum(num_mc1sf[1:3,:],axis=0)) / (den_mc1-np.sum(num_mc1[1:3,:],axis=0)) # pt dependent real_sf1
        real_sf1_err = np.sqrt(np.power(self.real_sf_err[0,:],2) + np.power(self.real_sf_err[1,:],2))
        self.real_sf = np.append([real_sf1],self.real_sf, axis=0)
        self.real_sf_err = np.append([real_sf1_err], self.real_sf_err, axis=0)
        num_mc1sf[0,:] = real_sf1*num_mc1[0,:]
        #
        num_mc, den_mc, eff_mc_up, eff_mc_down, edges = get_eff(self.mc2_df)
        eff_fake_mc = num_mc/den_mc
        #
        num_data, den_data, eff_data_up, eff_data_down, _ = get_eff(self.data_df)
        num_data_mc2 = num_data - num_mc1sf
        den_data_mc2 = den_data - den_mc1
        eff_fake_data = num_data_mc2/den_data_mc2
        #
        eff_fake_data_sys_up = abs(((num_data - 1.2*num_mc1sf)/(den_data - 1.2*den_mc1))-eff_fake_data)   # %20 norm unc
        eff_fake_data_sys_down = abs(((num_data - 0.8*num_mc1sf)/(den_data - 0.8*den_mc1))-eff_fake_data) # %20 norm unc
        #eff_fake_data_sys_sf_updown = self.real_sf_err*(den_mc1/den_mc)
        total_eff_fake_data_sys_up = np.sqrt(np.power(eff_fake_data_sys_up,2)+np.power(eff_data_up-(num_data/den_data),2))
        total_eff_fake_data_sys_down = np.sqrt(np.power(eff_fake_data_sys_down,2)+np.power(eff_data_down-(num_data/den_data),2))
        #print(eff_fake_mc)
        #
        fake_sf123 = eff_fake_data/eff_fake_mc
        fake_sf123_err_up   = fake_sf123*np.sqrt(np.power(total_eff_fake_data_sys_up/eff_fake_data,2) + np.power(abs(eff_mc_up- eff_fake_mc)/eff_fake_mc,2))
        fake_sf123_err_down = fake_sf123*np.sqrt(np.power(total_eff_fake_data_sys_down/eff_fake_data,2) + np.power(abs(eff_fake_mc - eff_mc_down)/eff_fake_mc,2))
        fake_sf123_err = np.max(np.vstack([[fake_sf123_err_up],[fake_sf123_err_down]]),axis=0)
        self.fake_sf = fake_sf123
        self.fake_sf_err = fake_sf123_err
        #

    def store_bbvl_sf(self):
        out_name = cfg.dataDir+f'/deepak8sf/deepak8_bbvL_sf_{self.year}.json'
        out_dict = {
            'score_bins':[0.]+self.score_bins,
            'pt_bins':self.pt_bins,
            'real': {
                'score_pt_sf'    : self.real_sf.tolist(),
                'score_pt_sf_err': self.real_sf_err.tolist(),
            },
            'fake': {
                'score_pt_sf'    : self.fake_sf.tolist(),
                'score_pt_sf_err': self.fake_sf_err.tolist(),
            },
                    
        }
        with open(out_name, 'w') as jsf:
            json.dump(out_dict, jsf, indent=4)
        
        
    def seperate_mc_dict(self):
        self.mc_dict   = {k: self.mc_dict[k][self.initial_datamc_cuts(self.mc_dict[k])] for k in self.mc_dict}
        nonhbb = self.mc_dict['ttH'].loc[self.mc_dict['ttH']['Hbb']==False,:].copy()
        hbb = self.mc_dict['ttH'].drop(nonhbb.index).copy()
        nonzbb = self.mc_dict['ttZ'].loc[self.mc_dict['ttZ']['Zbb']==False,:].copy()
        zbb = self.mc_dict['ttZ'].drop(nonzbb.index).copy()
        mc2_df = pd.concat([self.mc_dict[k] for k in ['TTBar','other','VJets','ttX','single_t']]+[nonzbb,nonhbb],axis='rows',ignore_index=True)
        mc1_df = pd.concat([self.mc_dict[k] for k in ['tt_bb','tt_2b']]+[zbb,hbb],axis='rows',ignore_index=True)
        del self.mc_dict
        mc1_df = self.add_sf_to_df(mc1_df)
        self.get_total_weight(mc1_df)
        self.get_total_weight(mc2_df)
        return mc1_df, mc2_df

    def add_sf_to_df(self, mc1_df):
        sf_dict = self.sf_dict['sig_sf']
        pt_bins, score_bins = sf_dict['pt_bins'], sf_dict['score_bins']
        self.pt_bins, self.score_bins = pt_bins, score_bins
        sf = sf_dict[self.year]['score_pt_sf']
        sf_err = sf_dict[self.year]['score_pt_sf_err']
        self.real_sf, self.real_sf_err = np.array(sf), np.array(sf_err)
        #
        pt_digi = pd.cut(mc1_df['Zh_pt'].clip(pt_bins[0]+0.001,pt_bins[-1]-0.001), bins = pt_bins, right=True,  include_lowest=True, labels=range(len(pt_bins)-1))
        score_digi = pd.cut(mc1_df['Zh_bbvLscore'], bins = score_bins, right=True,  include_lowest=True, labels=range(len(score_bins)-1))
        #pt_digi_ignore_nan = np.nan_to_num(pt_digi.values).astype(int)
        score_digi_ignore_nan = np.nan_to_num(score_digi.values).astype(int)
        #get_sf = (lambda _sf : np.array([ _sf[x][y] for x,y in zip(score_digi_ignore_nan,pt_digi_ignore_nan)]))
        get_sf = (lambda _sf : np.array([ _sf[x][y] for x,y in zip(score_digi_ignore_nan,pt_digi.values)]))
        mc1_df['dak8md_bbvl_sf'] = np.where(score_digi.isna().values, 1, get_sf(sf))
        mc1_df['dak8md_bbvl_sf_err'] = np.where(score_digi.isna().values, 0, get_sf(sf_err))
        #
        return mc1_df
        
        
    def get_sf_from_json(self):
        return json.load(open(self.bbvl_sf_file,'r'))

    def get_data_df(self):
        get_df = (lambda name: pd.read_pickle(f'{cfg.master_file_path}{self.year}/data_files/{name}_val.pkl'))
        data_df = pd.concat([get_df("Data_SingleElectron"),get_df("Data_SingleMuon")],axis='rows',ignore_index=True)
        #initial_data_cuts = (lambda df_:( ( (df_['isEleE']==True) | (df_['isMuonE']==True) ) &
        #                                  (df_['Zh_bbvLscore'] >= 0.0)))
        data_df = data_df[self.initial_datamc_cuts(data_df)]
        return data_df

    def get_total_weight(self, _df):
        #_df['tot_weight'] = 1
        _df['tot_weight'] = (_df['weight']* np.sign(_df['genWeight']) 
                             #* (np.where(self.mc_df['weight']>300,0,1))
                             * _df['topptWeight']
                             * (_df['HEM_weight']  if self.year == '2018' else 1.0 )
                             * _df['lep_trigeffsf']
                             * _df['lep_sf']
                             * _df['BTagWeight'] 
                             * _df['puWeight']  
                             * (_df['PrefireWeight'] if self.year != '2018' else 1.0))
        

if __name__ == '__main__':
    _ = DeepAk8_fakeSF('2016')
    _ = DeepAk8_fakeSF('2017')
    _ = DeepAk8_fakeSF('2018')
