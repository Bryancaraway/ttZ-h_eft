##### CALC LEPTON TRIG EFF  #####
### Written by: Bryan Caraway ###
#################################
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import uproot
import numpy as np
import awkward
import concurrent.futures
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
#import functools
import config.ana_cff as cfg
from modules.AnaDict import AnaDict
from modules.getdata import getData
from lib.fun_library import clop_pear_ci

class Calc_LepEff():
    '''
     calc lep eff per year
     this is done for pt bins
     and eta bins
    '''
    # coordinates to get data
    tag        = 'tight'
    roodir     = cfg.file_path
    treedir    = f'{cfg.tree_dir}_bb'
    data_tree = {'mu' :'MuData',
                 'ele':'EleData'}
    '''
    mu_id_map = {0:'loose',1:'med',2:'tight'}
    ele_id_map = {0:'fail', 1:'veto',2:'loose',3:'med',4:'tight'}
    '''
    # Object requirements for electron and muon
    iso_type = 'miniIso'

    mu_sel  = {'Id'     : 1,
               'Iso'    : 0.15,
               'miniIso': 0.2}

    ele_sel = {'Id'     : 4,
               'Iso'    : 9999,
               'miniIso': 0.1}

    lep_sel = {'mu':mu_sel, 'ele':ele_sel}
    # HLT path req per year and lepton
    hlt_path = cfg.hlt_path
    # Info to get from trees
    mu_keys   = ["Muon_pt","Muon_eta","Muon_phi","Muon_mass",
                 "Muon_miniPFRelIso_all","Muon_pfRelIso04_all",
                 "Muon_FlagId"]  # 0,1,2 = loose, med, tight
    ele_keys   = ["Electron_pt","Electron_eta","Electron_phi","Electron_mass",
                 "Electron_miniPFRelIso_all", 
                 "Electron_cutBasedNoIso", "Electron_cutBased"] # 0,1,2,3,4 = fail, veto, loose, med, tight
    hlt_keys  = cfg.ana_vars['dataHLT_all'] 
    # Pt, eta bins for lep eff caclulation
    pt_bins  = {'muon': [],
                'electron': []
    }
    eta_bins = {'muon': [-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,-0.2,
                         0.2,0.3,0.9,1.2,1.6,2.1,2.4],
                'electron': [-2.5,-1.566,-1.444,-0.8, 0.0, 
                             0.8,1.444,1.566,2.5]
    }
    #
    def __init__(self, year=None):
        self.year = year
        self.hlt_keys += cfg.ana_vars[f'dataHLT_{self.year}']
        self.roofile = f'Data_{self.year}_lep.root'
        #
        self.get_trees()
        self.calc_lep_eff('muon')
        self.calc_lep_eff('electron')
    
    def get_trees(self):
        f= uproot.open(self.roodir+self.roofile)
        self.t = {'mu' :  f.get(f"Training_bb/{self.data_tree['mu']}"),
                  'ele' : f.get(f"Training_bb/{self.data_tree['ele']}")}
        

    def calc_lep_eff(self,opt):
        '''
        Use SingleElectron[Muon] data to calculate muon[electron] trigger efficiency
        Need to pass analysis baseline with 
        at least one good electron AND one good muon
        store first good muon[electron] pt and eta 
        '''
        opt_dict = {'electron': self.t['mu'],
                    'muon':     self.t['ele']}
        tree = opt_dict[opt]
        #
        getData.DF_Container.set_attr(True, self.year, 3, 99, 1, None,None) # nominal, 3 jets with 1 bottom
        getData.DF_Container.set_current_tree_mask(tree)
        event_mask = getData.DF_Container.event_mask 
        #
        mu, ele, hlt = self.get_tree_data(tree)
        mu_iso = {'Iso'    : 'Muon_pfRelIso04_all',
                  'miniIso': 'Muon_miniPFRelIso_all'}
        ele_id  = {'Iso'    :'Electron_cutBased',
                  'miniIso':'Electron_cutBasedNoIso'}
        
        goodMu, goodEle, lep_event_mask = self.get_good_leps(mu, ele)
#            (mu[mu_iso[self.iso_type]]         < self.mu_sel[self.iso_type]),
#            (mu['Muon_FlagId']                >= self.mu_sel['Id']),
#            (ele['Electron_miniPFRelIso_all']  < self.ele_sel[self.iso_type]),
#            (ele[ele_id[self.iso_type]]       >= self.ele_sel['Id'])
#        )
        tot_event_mask = lep_event_mask & event_mask
        lep_pt_eta = {'electron': (ele['Electron_pt'][goodEle][tot_event_mask][:,0], 
                                   ele['Electron_eta'][goodEle][tot_event_mask][:,0]),
                      'muon'    : (mu['Muon_pt'][goodMu][tot_event_mask][:,0], 
                                   mu['Muon_eta'][goodMu][tot_event_mask][:,0])
        }
        lep_pt, lep_eta = lep_pt_eta[opt]
        hlt_mask = self.hlt_path[opt][self.year](hlt[tot_event_mask])
        # find best bins for pt and eta
        self.pt_bins[opt], _ = self.define_bins(lep_pt[hlt_mask], lep_eta[hlt_mask])
        #
        num_pt, num_pt_err, pt_edges    = self.histize(lep_pt[hlt_mask],  self.pt_bins[opt])
        num_eta, num_eta_err, eta_edges = self.histize(lep_eta[hlt_mask], self.eta_bins[opt])
        #
        den_pt, den_pt_err,  _  = self.histize(lep_pt,  self.pt_bins[opt])
        den_eta, den_eta_err, _ = self.histize(lep_eta, self.eta_bins[opt])
        


        eff_pt     = num_pt/den_pt
        #eff_pt_err = abs(eff_pt)*np.sqrt( np.power( num_pt_err/num_pt ,2) + np.power( den_pt_err/den_pt ,2))
        eff_pt_err = clop_pear_ci(num_pt,den_pt,return_error=True)
        eff_pt_down, eff_pt_up = clop_pear_ci(num_pt,den_pt)
        
        eff_eta = num_eta/den_eta
        #eff_eta_err = abs(eff_eta)*np.sqrt( np.power( num_eta_err/num_eta ,2) + np.power( den_eta_err/den_eta ,2))
        eff_eta_err = clop_pear_ci(num_eta,den_eta,return_error=True)
        eff_eta_down, eff_eta_up = clop_pear_ci(num_eta,den_eta)
        #
        print(num_pt,'\n',den_pt,'\n')
        print(num_eta,'\n',den_eta,'\n')
        print(self.pt_bins[opt])
        print(eff_pt,'\n' ,eff_pt_err, '\n')
        #self.plot_num_den(num_pt, num_pt_err, den_pt, den_pt_err, pt_edges, opt)
        #self.plot_eff(eff_pt, eff_pt_err, pt_edges, eff_eta, eff_eta_err, eta_edges, opt)
        self.eff_to_pickle(eff_pt, eff_pt_up, eff_pt_down, eff_eta, eff_eta_up, eff_eta_down, opt)
        
    def get_tree_data(self,tree):
        executor = concurrent.futures.ThreadPoolExecutor()
        #
        mu_info   = AnaDict({k:tree.array(k, executor=executor) for k in self.mu_keys})
        ele_info  = AnaDict({k:tree.array(k, executor=executor) for k in self.ele_keys})
        hlt_info  = AnaDict({k:tree.array(k, executor=executor) for k in self.hlt_keys})
        #
        #mu_info  = mu_info[ (mu_info['Muon_pt']      >= 30) & (abs(mu_info['Muon_eta'])      <= 2.4)]
        #ele_info = ele_info[(ele_info['Electron_pt'] >= 30) & (abs(ele_info['Electron_eta']) <= 2.5)]
        #
        return mu_info, ele_info, hlt_info
    
    def eff_to_pickle(self, eff_pt, eff_pt_up, eff_pt_down, eff_eta, eff_eta_up, eff_eta_down, lep_type):
        eff = AnaDict({'pt':{}, 'eta':{}})
        for i in range(len(eff_pt)):
            eff['pt'][f'{self.pt_bins[lep_type][i]},{self.pt_bins[lep_type][i+1]}'] = {
                'values':eff_pt[i].round(3), 
                'up':eff_pt_up[i].round(3), 
                'down':eff_pt_down[i].round(3)}
        #
        for i in range(len(eff_eta)):
            eff['eta'][f'{self.eta_bins[lep_type][i]},{self.eta_bins[lep_type][i+1]}'] = {
                'values':eff_eta[i].round(3), 
                'up':eff_eta_up[i].round(3), 
                'down':eff_eta_down[i].round(3)}

        eff_dir  = f'files/{self.year}/eff_files/'
        pkl_file = f'{lep_type}_eff_{self.year}_{self.tag}.pkl'
        eff.to_pickle(eff_dir+pkl_file)
        

    def plot_eff(self, eff_pt, err_pt, edges_pt, eff_eta, err_eta, edges_eta, lep_type):
        fig, ax = plt.subplots(1,2, figsize=(16,10))
        for i, (eff, err, edges) in enumerate([(eff_pt,err_pt,edges_pt),(eff_eta,err_eta,edges_eta)]):
            
            bin_c = (edges[1:]+edges[:-1])/2
            bin_w = (edges[1:]-edges[:-1])/2

            ax[i].errorbar(x=bin_c, xerr=bin_w,
                           y=eff,   yerr=err,
                           fmt='.', barsabove=True, capsize=2)
            ax[i].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i].yaxis.set_minor_locator(AutoMinorLocator())
            ax[i].set_ylim(0, 1.1)
            ax[i].grid(True)
        plt.grid(True)
        fig.suptitle(f'{lep_type} trigger eff. for Data, {self.year}')
        plt.show()
        plt.clf()
        plt.close()

    def plot_num_den(self, num, num_err, den, den_err, edges, lep_type):
        fig, ax = plt.subplots(1,1, figsize=(16,10))

        bin_c = (edges[1:]+edges[:-1])/2
        bin_w = (edges[1:]-edges[:-1])/2

        ax.errorbar(x=bin_c, xerr=bin_w,
                    y=num,   yerr=num_err,
                    fmt='.', barsabove=True, capsize=2,
                    label='num')
        ax.errorbar(x=bin_c, xerr=bin_w,
                    y=den,   yerr=den_err,
                    fmt='.', barsabove=True, capsize=2,
                    label='den')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True)
        ax.legend()
        fig.suptitle(f'{lep_type} yield for Data, {self.year}')
        plt.show()
        plt.clf()
        plt.close()

    @staticmethod
    def histize(x, bins):
        n_, edges_ = np.histogram(x, bins=bins, range=(bins[0],bins[-1]))
        return n_, np.sqrt(n_), edges_


    def get_good_leps(self, mu_info, el_info):#mu_iso_cut, mu_id_cut, ele_iso_cut, ele_id_cut):
        #mu_mask  = ((mu_iso_cut)  & (mu_id_cut))
        #ele_mask = ((ele_iso_cut) & (ele_id_cut))
        mu_mask = cfg.lep_sel['muon'](mu_info)
        ele_mask = cfg.lep_sel['electron'][self.year](el_info)
        #
        lep_event_mask = ((mu_mask[mu_mask].counts > 0) & (ele_mask[ele_mask].counts > 0))
        return mu_mask, ele_mask, lep_event_mask
    
    @staticmethod
    def define_bins(pt_, eta_):
        quants = np.linspace(0.0, 1.0, num=6)
        #print(quants)
        pt_bins = np.append(np.quantile(pt_[pt_< 200], quants).round(0), [300,700])
        eta_bins = [-2.4,-1.6, -0.8, 0.0, 0.8, 1.6, 2.4]
        return pt_bins , eta_bins

    @classmethod
    def set_lep_sel(cls, lep_sel_dict):
        cls.lep_sel = lep_sel_dict
        cls.mu_sel  = lep_sel_dict['mu']
        cls.ele_sel = lep_sel_dict['ele']

if __name__ == '__main__':
    #
    eff = Calc_LepEff(year=sys.argv[1])
