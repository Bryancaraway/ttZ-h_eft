import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import uproot
import concurrent.futures
executor = concurrent.futures.ThreadPoolExecutor()
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from lib.fun_library import fill1e, fillne, t2Run, save_pdf, deltaR
from modules.AnaDict import AnaDict
import config.ana_cff as cfg

roo_dir = 'roc_roots/'
roo_files = ['tt_sl.root','tthbb.root','ttbb_sl.root','ttzqq.root']

class Test_NanoAODv7:

    gen_vars = cfg.ana_vars['genpvars']
    fj_vars  = [
        'FatJet_msoftdrop','FatJet_pt','FatJet_eta','FatJet_phi','FatJet_btagHbb',
        'FatJet_btagDDBvL','FatJet_deepTagMD_ZHbbvsQCD','FatJet_deepTagMD_bbvsLight']
    
    def __init__(self,roo_files=roo_files, roo_dir=roo_dir): # initialize files which to compare
        self.roo_files = roo_files
        self.roo_dir   = roo_dir
        self.data = {}
        self.df = pd.DataFrame()
        self.getdata()
        #
        self.process_sig('tthbb')
        self.process_sig('ttzqq')
        #
        self.process_bkg('tt_sl')
        self.process_bkg('ttbb_sl')
        #
        self.apply_weights()
        self.normalize_scores()
        #
        #print(self.df[self.df['is_tt']==True])
        #print(self.df[self.df['is_tt_bb']==True])
        #print(self.df[self.df['matchedGen_Zbb']==True])
        #print(self.df[self.df['matchedGen_Hbb']==True])
        #
        self.mvas = ['best_fj_hbb','best_fj_dd_bvl','best_fj_dt_zhbb_md','best_fj_dt_bbvl_md']
        #plt.hist(self.df['best_fj_dd_bvl'][self.df['is_tt'] == True], density=True)
        #plt.hist(self.df['best_fj_dt_bbvl_md'][self.df['Hbb'] == True], density=True)
        #plt.show()
        #print(self.df)

    @t2Run
    def getdata(self):
        for f in self.roo_files:
            s = f.replace('.root','')
            self.data[s] = {}
            t = uproot.open(roo_dir+f)['Events']
            self.data[s]['gen']     = AnaDict({v: t.array(v, executor=executor) for v in self.gen_vars})
            self.data[s]['fj'] = AnaDict({v: t.array(v, executor=executor) for v in self.fj_vars}) 
            # exlude fj outside of zh reconstruction potential
            fj_cut = ((self.data[s]['fj']['FatJet_pt'] > 300) & (abs(self.data[s]['fj']['FatJet_eta']) <= 2.4) & 
                   (self.data[s]['fj']['FatJet_msoftdrop'] >= 50) & (self.data[s]['fj']['FatJet_msoftdrop'] <= 200))
            for fj_v in self.data[s]['fj']:
                self.data[s]['fj'][fj_v] = self.data[s]['fj'][fj_v][fj_cut][fj_cut.sum() > 0]
            for gen_v in self.data[s]['gen']:
                self.data[s]['gen'][gen_v] = self.data[s]['gen'][gen_v][fj_cut.sum() > 0]
            
    @t2Run
    def process_bkg(self,b):
        df = pd.DataFrame()
        gen_tt_bb = self.data[b]['gen']['genTtbarId'] % 100
        is_tt_bb = gen_tt_bb >= 51
        df['is_tt'] = ((~is_tt_bb) & (b == 'tt_sl'))
        df['is_tt_bb'] = ((is_tt_bb) & (b == 'ttbb_sl'))
        # need to get the best scoring candidate in the event
        fj_pt, fj_eta, fj_phi, fj_sdm = map(fillne, self.data[b]['fj'].loc(['FatJet_pt','FatJet_eta','FatJet_phi', 'FatJet_msoftdrop']).values())
        fj_kinem_keys = ['pt','eta','phi','sdm']
        fj_hbb, fj_dd_bvl, fj_dt_zhbb_md, fj_dt_bbvl_md = map(fillne, self.data[b]['fj'].loc(['FatJet_btagHbb','FatJet_btagDDBvL','FatJet_deepTagMD_ZHbbvsQCD','FatJet_deepTagMD_bbvsLight']).values())
        fj_mva_keys = ['best_fj_hbb','best_fj_dd_bvl','best_fj_dt_zhbb_md','best_fj_dt_bbvl_md']
        #
        #
        for mva, k in zip([fj_hbb, fj_dd_bvl, fj_dt_zhbb_md, fj_dt_bbvl_md], fj_mva_keys):
            best_ind  = -1*np.argsort(-1*mva, axis=1) # have to do -1*(-1) to get max in first index
            get_best = (lambda fj_v: np.take_along_axis(fj_v, best_ind, axis=1)[:,0])
            df[k] = get_best(mva)
            #for kinem, kinem_key in zip([fj_pt, fj_eta, fj_phi, fj_sdm], fj_kinem_keys):
            #    pass
                #df[f'{k}_{kinem_key}'] = get_best(kinem)
            #
        #
        self.df = pd.concat([self.df,df.dropna(axis='rows')], ignore_index=True)

    
    @t2Run
    def process_sig(self,s):
        df = pd.DataFrame()
        g = 'GenPart_'
        gen_ids, gen_mom, gen_pt, gen_eta, gen_phi, gen_mass = self.data[s]['gen'].loc(
            [g+'pdgId',g+'genPartIdxMother',g+'pt',g+'eta',g+'phi',g+'mass']
        ).values()
        islep   = (((abs(gen_ids) == 11) | (abs(gen_ids) == 13)) 
                   & ((abs(gen_ids[gen_mom[gen_mom]]) == 6) & (abs(gen_ids[gen_mom]) ==24)))
        #
        fj_pt, fj_eta, fj_phi, fj_sdm = map(fillne, self.data[s]['fj'].loc(['FatJet_pt','FatJet_eta','FatJet_phi', 'FatJet_msoftdrop']).values())
        #
        isbb_fromZ     = ((abs(gen_ids) == 5) & (gen_ids[gen_mom] == 23))
        isbb_fromH    = ((abs(gen_ids) == 5) & (gen_ids[gen_mom] == 25))
        isHbb     = ((gen_ids == 25) & (isbb_fromH.sum() == 2))
        isZbb  = ((gen_ids == 23) & (isbb_fromZ.sum() == 2) & ((isHbb.sum() == 0)))
        isZH = ((isHbb) | (isZbb))
        #
        fill_cut = (lambda g: fill1e(g[isZH]).flatten())
        zh_pt, zh_eta, zh_phi = map(fill_cut, [gen_pt, gen_eta, gen_phi])
        #
        #zh_match_dR = deltaR(zh_eta,zh_phi,rZh_eta, rZh_phi)
        zh_match_dR = deltaR(zh_eta,zh_phi,fj_eta,fj_phi)
        ind_genm_fj = np.argsort(zh_match_dR,axis=1)
        first_closest = (lambda fj_mva: np.take_along_axis(fillne(fj_mva),ind_genm_fj,axis=1)[:,0])
        fj_hbb, fj_dd_bvl, fj_dt_zhbb_md, fj_dt_bbvl_md = map(
            first_closest, 
            self.data[s]['fj'].loc(['FatJet_btagHbb','FatJet_btagDDBvL','FatJet_deepTagMD_ZHbbvsQCD','FatJet_deepTagMD_bbvsLight']).values())
        nfj_pt, nfj_eta, nfj_phi, nfj_sdm, nfj_genm_dr = map((lambda fj_v: np.take_along_axis(fj_v,ind_genm_fj,axis=1)[:,0]), [fj_pt,fj_eta,fj_phi,fj_sdm, zh_match_dR])
        #df['best_fj_hbb_pt'] = df['best_fj_dd_bvl_pt'] = df['best_fj_dt_zhbb_md_pt'] = df['best_fj_dt_bbvl_md_pt'] = nfj_pt
        #df['best_fj_hbb_eta'] = df['best_fj_dd_bvl_eta'] = df['best_fj_dt_zhbb_md_eta'] = df['best_fj_dt_bbvl_md_eta'] = nfj_eta
        #df['best_fj_hbb_phi'] = df['best_fj_dd_bvl_phi'] = df['best_fj_dt_zhbb_md_phi'] = df['best_fj_dt_bbvl_md_phi'] = nfj_phi
        #df['best_fj_hbb_sdm'] = df['best_fj_dd_bvl_sdm'] = df['best_fj_dt_zhbb_md_sdm'] = df['best_fj_dt_bbvl_md_sdm'] = nfj_sdm
        #df['fj_genm_dr'] = nfj_genm_dr
        #
        df['best_fj_hbb'] = fj_hbb
        df['best_fj_dd_bvl'] = fj_dd_bvl
        df['best_fj_dt_zhbb_md'] = fj_dt_zhbb_md
        df['best_fj_dt_bbvl_md'] = fj_dt_bbvl_md
        #
        zh_match = ( (zh_pt >= (150.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
        zh_match = zh_match[((fj_pt[zh_match_dR <= 0.9] >= 200 ).sum() >= 1)]
        #zh_match = ((zh_match_dR <= 0.8) & (zh_pt >= ( self.pt_cut-100.)) & (zh_eta <= 2.4) & (zh_eta >= -2.4))
        #
        df['Zbb']= (isZbb.sum() > 0)
        df['Hbb']= (isHbb.sum() > 0)
        good_lep = ( (gen_pt[islep] > 30) & (abs(gen_eta[islep]) < 2.5) )
        df['matchedGenLep']   = (good_lep.sum() == 1)
        df['matchedGen_Zbb']  = (((zh_match).sum() > 0) & (df['matchedGenLep']) & (isZbb.sum() >  0))
        df['matchedGen_Hbb']  = (((zh_match).sum() > 0) & (df['matchedGenLep']) & (isHbb.sum() >  0))
        #
        self.df = pd.concat([self.df,df.dropna(axis='rows')], ignore_index=True)

    def apply_weights(self):
        # define signal and bkg
        self.df['Y'] = 1
        self.df['Y'].loc[self.df['is_tt_bb']==True] = 0
        self.df['Y'].loc[self.df['is_tt']==True] = 0
        self.df.drop(index=self.df[(self.df['Hbb'] == True) & (self.df['matchedGen_Hbb'] == False)].index , inplace=True)
        self.df.drop(index=self.df[(self.df['Zbb'] == True) & (self.df['matchedGen_Zbb'] == False)].index , inplace=True)
        #
        w_dict = {'ttbb_sl': (3771+1604),
                  'tt_sl'  : 18269.9,
                  'tthbb'  : 195,
                  'ttzbb'  : 80,
        }
        self.df['weight'] = 1
        self.df['weight'].loc[self.df['is_tt_bb']==True] = w_dict['ttbb_sl']/len(self.df['is_tt_bb']==True)
        self.df['weight'].loc[self.df['is_tt']==True]    = w_dict['tt_sl']  /len(self.df['is_tt']==True)
        self.df['weight'].loc[self.df['Zbb']==True]      = w_dict['ttzbb']  /len(self.df['Zbb']==True)
        self.df['weight'].loc[self.df['Hbb']==True]      = w_dict['tthbb']  /len(self.df['Hbb']==True)
        #normalize to highest weight
        self.df.loc[:,'weight'] = self.df['weight']/min(self.df['weight'])
        
    def normalize_scores(self):
        
        for mva in ['best_fj_hbb','best_fj_dd_bvl','best_fj_dt_zhbb_md','best_fj_dt_bbvl_md']:
            df = self.df[mva][self.df[mva] > -1]
            self.df.loc[self.df[mva] > -1 , mva] = (df - min(df)) / (max(df) - min(df))

    @save_pdf('test_zhbb_mvas.pdf')
    def make_rocplots(self):
        # first compare against ttbar
        sb = self.df[(self.df['is_tt'] == True) | (self.df['Hbb'] == True) | (self.df['Zbb'] == True)]
        self.plot_roc(sb, 'MVA comparisons (sig vs tt)')
        sb = self.df[(self.df['is_tt_bb'] == True) | (self.df['Hbb'] == True) | (self.df['Zbb'] == True)]
        self.plot_roc(sb, 'MVA comparisons (sig vs ttbb)')

    def plot_roc(self,df,title):
        from sklearn import metrics
        fig, ax = plt.subplots()
        ax.plot([0,1],[0,1], 'k--')
        #cut = (lambda m: (self.df[f'{m}_pt'] > 200) & (abs(self.df[f'{m}_eta']) <= 2.4) & 
        #       (self.df[f'{m}_sdm'] >= 50) & (self.df[f'{m}_sdm'] <= 200) & (self.df[m] > -1))
        roc_dict = {}
        for mva in self.mvas:
            df = df[df[mva] > -1]
            #print(df[mva],'\n',min(df[mva]),max(df[mva]))
            fpr, tpr, thresh = metrics.roc_curve(df['Y'], df[mva], sample_weight=df['weight'])
            roc_dict[mva] = {'tpr':tpr,'fpr':fpr}
            ax.plot(fpr,tpr,
                    label=f'{mva}, AUC: {metrics.auc(fpr,tpr):.4f}')
        #
        if len(roc_dict['best_fj_dt_bbvl_md']['fpr']) >= len(roc_dict['best_fj_hbb']['fpr']):
            fpr_coor = roc_dict['best_fj_dt_bbvl_md']['fpr']
            new_tpr = np.interp(roc_dict['best_fj_dt_bbvl_md']['fpr'], roc_dict['best_fj_hbb']['fpr'], roc_dict['best_fj_hbb']['tpr'])
            diff_tpr = roc_dict['best_fj_dt_bbvl_md']['tpr'] - new_tpr
            div_tpr  = roc_dict['best_fj_dt_bbvl_md']['tpr']/new_tpr
        else:
            fpr_coor = roc_dict['best_fj_hbb']['fpr']
            new_tpr = np.interp(roc_dict['best_fj_hbb']['fpr'], roc_dict['best_fj_dt_bbvl_md']['fpr'], roc_dict['best_fj_dt_bbvl_md']['tpr'])
            diff_tpr = new_tpr - roc_dict['best_fj_hbb']['tpr']
            div_tpr = new_tpr / roc_dict['best_fj_hbb']['tpr']
        #
        ax.plot(
            fpr_coor,
            diff_tpr,
            label = 'dt_bbvl_md tpr - best_fj_hbb tpr (left axis)',
            color='k',
        )
        ax2 = ax.twinx()
        ax2.plot(fpr_coor,
                 div_tpr,
                 label=' dt_bbvl_md tpr/ best_fj_hbb tpr (right axis)',
                 color = 'gold')
        ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax2.legend(fontsize='x-small')
        #
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_ylim(0,1.05)
        ax.set_xlim(0,1)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, which='both')
        fig.suptitle(title)
        ax.legend(fontsize='x-small')
        
        #plt.show()


def main():
    _ = Test_NanoAODv7()
    #save_pdf('test_zhbb_mvas.pdf')(_.plot_roc)('MVA comparisons')
    _.make_rocplots()

if __name__ == '__main__' :
    main()
