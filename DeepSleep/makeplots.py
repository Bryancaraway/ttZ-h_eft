from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import numpy as np
from modules.AnaDict import AnaDict
from lib.fun_library import save_pdf, getFakebbvlCuts, getFakebbvlWeights
import config.ana_cff as cfg
###----=======-----###

#sepGenOpt ; 'sepGenSig','sepGenBkg','sepGenMatchedSig','sepGenMatchedBkg'

###----=======-----###

#jjec = 'ak8'
#jec_list = ['JESUp','JESDown','JERUp','JERDown']
#processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets','other']
#processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets']
processes = ['ttZ','ttH','TTBar','tt_B','ttX','single_t','VJets']
#processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b']

#@save_pdf('control_plots_anastrat.pdf')
#@save_pdf('control_plots_genzhpt.pdf')
#@save_pdf('ttzh_purity.pdf')
#@save_pdf('Zh_info.pdf')

#@save_pdf('control_plots_tight.pdf')


@save_pdf('ttbb_lheht.pdf')


#@save_pdf('NN_compare.pdf')
#@save_pdf('fakebbvlsf_CR.pdf')
#@save_pdf('fakebbvlsf_CR_Andrew.pdf')

#@save_pdf("rare_yields.pdf")
#@save_pdf("met_withqcd.pdf")
##@save_pdf("nlepfrac_study.pdf")

#@save_pdf('btv_study.pdf')

#@save_pdf('mass_sensitivity.pdf')


def main():
    for y in cfg.Years: 
        #for y in ['2018']: 
        print(y)
        #for jec in jec_list:
        #Plotter.load_data(y, addBSF=False, tag=f'{jjec}{jec}') #tag='ak4JESUp'
        #Plotter.load_data(y, samples=cfg.Sig_MC+cfg.Bkg_MC, addBSF=False, byprocess=True)
        Plotter.load_data(y, samples=['ttbb'], addBSF=False, byprocess=True)
        #Plotter.load_data(y, samples=cfg.Sig_MC+cfg.Bkg_MC+["QCD"], addBSF=False, byprocess=True)
        ''' LOOK AT STACKED DATA VS MC '''

        # --- control plots tight
        #StackedHist(processes,    'PV_npvsGood', xlabel= r'nPVs', bin_range=[0,70],    n_bins=35,   doCuts=True,  addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'Lep_pt', xlabel= r'lepton $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,   doCuts=True,  addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'lepton $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_pt', xlabel= r'Electron $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepElec==1', doCuts=True, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'Electron $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepElec==1',  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_pt', xlabel= r'Muon $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepMu==1', doCuts=True, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'Muon $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepMu==1',  doCuts=True, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'MET_pt', xlabel=r'missing $e_{T}$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'MET_phi', xlabel=r'missing $\phi$ (GeV)', bin_range=[-3.15,3.15],  n_bins=20,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'nBottoms', xlabel='# of ak4 b-jets', bin_range=[-0.5,8.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],  doCuts=True, addData=True, doShow=False) 
        #StackedHist(processes,    'n_ak4jets', xlabel='# of ak4 jets',  bin_range=[-0.5,12.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5],  doCuts=True, addData=True, doShow=False) 
        #StackedHist(processes,    'n_ak8jets', xlabel='# of ak8 jets', bin_range=[-0.5,6.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5],  doCuts=True, addData=True, doShow=False) 
        ##
        #StackedHist(processes,    'jetpt_1', xlabel=r'leading jet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'jetpt_2', xlabel=r'next-to-leading jet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetpt_1', xlabel=r'leading bjet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetpt_2', xlabel=r'next-to-leading bjet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'fjetpt_1', xlabel=r'leading fatjet $p_{T} (AK8)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'jeteta_1', xlabel=r'leading jet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'jeteta_2', xlabel=r'next-to-leading jet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'bjeteta_1', xlabel=r'leading bjet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'bjeteta_2', xlabel=r'next-to-leading bjet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'fjeteta_1', xlabel=r'leading fatjet $\eta (AK8)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'jetbtag_1', xlabel=r'leading jet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'jetbtag_2', xlabel=r'next-to-leading jet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetbtag_1', xlabel=r'leading bjet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetbtag_2', xlabel=r'next-to-leading bjet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'fjetsdm_1', xlabel=r'leading fatjet $m_{sd}$ (GeV)', bin_range=[50,200],  n_bins=25,     doCuts=True, addData=True, doShow=False)  

        ###########
        #StackedHist(processes,    'fjetbbvl_1', xlabel=r'leading fatjet deepAK8MD bbvL', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        
        ###   #
        #StackedHist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='Zh_bbvLscore>0.8',  doCuts=True, addData=True, doShow=False)  
        #Hist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ [GeV] (Zhpt>300,NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='NN>0.8;Zh_pt>300',  doCuts=True, addSoB=True, doShow=False)  
        #Hist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ [GeV] (Zhpt>450,NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='NN>0.8;Zh_pt>450',  doCuts=True, addSoB=True, doShow=False)  
        ##
        ##StackedHist(processes,    'Lep_pt', xlabel= r'lepton $p_{T}$ (GeV) lep_eta>2, tight baseline', bin_range=[20,750],    n_bins=30, add_cuts='Lep_eta>2;MET_pt>20',   doCuts=True, addData=True,  doShow=False)  
        ##StackedHist(processes,    'Lep_eta', xlabel=r'lepton $\eta$, tight baseline',        bin_range=[-2.6,2.6],  n_bins=26,   doCuts=True, add_cuts='MET_pt>20', addData=True, doShow=False)  
        ##StackedHist(processes,    'Lep_pt', xlabel= r'Electron $p_{T}$ (GeV) lep_eta>2. tight baseline', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepElec==1;Lep_eta>2;MET_pt>20', doCuts=True, addData=True,  doShow=False)  
        ##StackedHist(processes,    'Lep_eta', xlabel=r'Electron $\eta$ tight baseline',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepElec==1;MET_pt>20',  doCuts=True, addData=True, doShow=False)  
        ##StackedHist(processes,    'Lep_pt', xlabel= r'Muon $p_{T}$ (GeV) lep_eta>2, tight baseline', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepMu==1;Lep_eta>2;MET_pt>20', doCuts=True, addData=True,  doShow=False)  
        ##StackedHist(processes,    'Lep_eta', xlabel=r'Muon $\eta$ tight baseline',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepMu==1;MET_pt>20',  doCuts=True, addData=True, doShow=False)  
        ## ----
        ##StackedHist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLscore (no extra cuts)', bin_range=[0,1],  n_bins=20,  doCuts=False, addData=True, doShow=False)  
        ##StackedHist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLscore Fake?', bin_range=[0,1],  n_bins=20, add_cuts='Zh_M>65;Zh_M<90;Zh_closeb_invM>150;Zh_closeb_invM<190',  doCuts=False, addData=True, doShow=False)  
        ##StackedHist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLscore Fake? (>2 b-jets outZH)', bin_range=[0,1],  n_bins=20, add_cuts='Zh_M>65;Zh_M<90;Zh_closeb_invM>150;Zh_closeb_invM<190;n_b_outZh>=2',  doCuts=False, addData=True, doShow=False)  
        #
        #StackedHist(processes,'withdak8md_NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='withdak8md_NN<=0.7',  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'noak8md_NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='noak8md_NN<=0.7',     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'withbbvl_NN', bin_range=[0,1],  n_bins=20,  add_d_cuts='withbbvl_NN<=0.7',     doCuts=True, addData=True, doShow=False)  

        #StackedHist(processes,'NN', xlabel='new_withbbvl_NN (Zhpt1)', bin_range=[0,1], bins=[0. ,        0.47342254, 0.7280452,  0.77745633, 0.83437309, 0.89884436, 1],  n_bins=20, add_cuts='Zh_pt>200;Zh_pt<300',      doCuts=True, addData=False, doShow=True) 
        #StackedHist(processes,'NN', xlabel='new_withbbvl_NN (Zhpt2)', bin_range=[0,1], bins=[0. ,        0.47342254, 0.7280452,  0.77745633, 0.83437309, 0.89884436, 1],  n_bins=20, add_cuts='Zh_pt>300;Zh_pt<450',      doCuts=True, addData=False, doShow=False) 
        #StackedHist(processes,'NN', xlabel='new_withbbvl_NN (Zhpt3)', bin_range=[0,1], bins=[0. ,        0.47342254, 0.7280452,  0.77745633, 0.83437309, 0.89884436, 1],  n_bins=20, add_cuts='Zh_pt>450',                doCuts=True, addData=False, doShow=False) 
        #
        ## --- Study for Andrew
        #StackedHist(processes,'NN', xlabel='NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='NN<=0.4',     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'NN', xlabel='NN (Zh_pt>300)', bin_range=[0,1],  n_bins=20, add_cuts='Zh_pt>300', add_d_cuts='NN<=0.4',     doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'NN', xlabel='NN (Zh_pt>300, ZH M bins)', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80;Zh_pt>300;Zh_M>80;Zh_M<155', add_d_cuts='NN<=0.4',     doCuts=True, addData=True, doShow=False)  
        #Hist(['ttZ','ttH'],'NN', xlabel='NN (Zh_pt>300, ZH M bins)', bin_range=[0,1],  n_bins=20, add_cuts='NN>=.80;Zh_pt>300;Zh_M>80;Zh_M<155', doNorm=False, doCuts=True, doShow=False)  
        #Hist(['ttX'],'NN',        xlabel='NN (Zh_pt>300, ZH M bins)', bin_range=[0,1],  n_bins=20, sepGenOpt='sepBySample', add_cuts='NN>=.80;Zh_pt>300;Zh_M>80;Zh_M<155', doNorm=False,  doCuts=True, doShow=False)  
        #
        ### QCD study
        #Hist(['QCD','ttZ','ttH'],'min_sel_soft_mu_invm',xlabel='min_sel_soft_mu_invm', bin_range=[0,100],  n_bins=100, sepGenOpt='sepGenMatchedSig;+',   add_cuts='passSingleLepMu==1;NN<=1.80',    doCuts=True, doLog=True, doNorm=False, doShow=False)
        ###
        #Hist(['QCD','ttZ','ttH'],'Muon_sip3d',xlabel='Muon_sip3d', bin_range=[0,10],  n_bins=40, sepGenOpt='sepGenMatchedSig;+',  add_cuts='passSingleLepMu==1;NN<=1.80;passNotHadLep==1',    doCuts=True, doLog=True, doNorm=False, doShow=False)
        ###
        #
        #Hist(['QCD','ttZ','ttH'],'min_sel_soft_elec_invm',xlabel='min_sel_soft_elec_invm', bin_range=[0,100],  n_bins=100, sepGenOpt='sepGenMatchedSig;+',   add_cuts='passSingleLepElec==1;NN<=1.80',    doCuts=True, doLog=True, doNorm=False, doShow=False)
        #Hist(['QCD','ttZ','ttH'],'Electron_sip3d',xlabel='Electron_sip3d', bin_range=[0,10],  n_bins=40, sepGenOpt='sepGenMatchedSig;+',  add_cuts='passSingleLepElec==1;NN<=1.80;Zh_pt>300;passNotHadLep==1',    doCuts=True, doLog=True, doNorm=False, doShow=False)
        ###
        #
        ##
        #StackedHist(processes,'min_sel_soft_elec_invm', bin_range=[0,100], n_bins=100, add_cuts='passSingleLepElec==1', doCuts=True, doLog=True, addData=True, doShow=False)
        #StackedHist(processes,'min_sel_soft_mu_invm', xlabel='min_sel_soft_mu_invm', bin_range=[0,100], n_bins=100, add_cuts='passSingleLepMu==1', doCuts=True, doLog=True, addData=True, doShow=False)
        #Hist(['QCD','ttZ','ttH'],'MET_pt', xlabel=r'missing $e_{T}$ (GeV), with j/psi lepton cut', bin_range=[0,500],  n_bins=25,  add_cuts='passNotHadLep==1', sepGenOpt='sepGenMatchedSig;+',   doCuts=True, doShow=False)  
        #Hist(['QCD','ttZ','ttH'],'MET_pt', xlabel=r'missing $e_{T}$ (GeV), with j/psi lepton cut, Zhpt>300', bin_range=[0,500],  n_bins=25,  add_cuts='passNotHadLep==1;Zh_pt>300', sepGenOpt='sepGenMatchedSig;+',   doCuts=True, doShow=False)  
        #StackedHist(processes+['QCD'],'MET_pt', xlabel=r'missing $e_{T}$ (GeV), with j/psi lepton cut and sip3d < 4', bin_range=[0,500],  n_bins=25,  doCuts=True, addData=True, doShow=False)  

        # ------------- N lep fraction -------- #

        #Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps', xlabel='n_tt_leps', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, doNorm=True, doShow=False)
        #Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps', xlabel='n_tt_leps (NN>0.8)', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, add_cuts='NN>0.8', doNorm=True, doShow=False)
        #
        #Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps_notau', xlabel='n_tt_leps no taus', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, doNorm=True, doShow=False)
        #Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps_notau', xlabel='n_tt_leps no taus (NN>0.8)', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, add_cuts='NN>0.8', doNorm=True, doShow=False)

        # --------------------------- #
        

        #Hist(['ttZ','ttH'],'NN',xlabel='NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', sepGenOpt='sepGenSig', doCuts=True, doLog=True, doNorm=False, doShow=True)  
        #Hist(['ttZ','ttH'],'withbbvl_NN',xlabel='NN (Zh_pt>300)', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80;Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH','tt_B','TTBar'],'NN',xlabel=r'DNN score (${p}_{\mathrm{T}}^{\mathrm{Z/H}}>300$ GeV)', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80;Zh_pt>300', sepGenOpt='sepGenMatchedSig;+', doCuts=True, doLog=True, doNorm=True, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_bbvLscore',xlabel='bbvL', bin_range=[.8,1],  n_bins=5, add_cuts='NN<=1.80', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_bbvLscore',xlabel='bbvL (Zh_pt>300)', bin_range=[.8,1],  n_bins=5, add_cuts='NN<=1.80;Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_bbvLscore',xlabel='bbvL (Zh_pt>450)', bin_range=[.8,1],  n_bins=5, add_cuts='NN<=1.80;Zh_pt>450', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_M',xlabel='Z/H $m_{sd}$ ',   bins=[50,75,90,105,120,140,200],  sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_M',xlabel='Z/H $m_{sd}$ (Zh_pt>300)', bins=[50,75,90,105,120,140,200], add_cuts='Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_M',xlabel='Z/H $m_{sd}$ (Zh_pt>450)', bins=[50,75,90,105,120,140,200], add_cuts='Zh_pt>450', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_closeb_dr',xlabel='$\Delta R$(Z/H, closest b) ',          bin_range=[0,3], bins=15,  sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_closeb_dr',xlabel='$\Delta R$(Z/H, closest b) (Zh_pt>300)', bin_range=[0,3], bins=15, add_cuts='Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_closeb_dr',xlabel='$\Delta R$(Z/H, closest b) (Zh_pt>450)', bin_range=[0,3], bins=15, add_cuts='Zh_pt>450', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_2ndcloseb_dr',xlabel='$\Delta R$(Z/H, 2nd closest b) ',          bin_range=[0,3], bins=15,  sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_2ndcloseb_dr',xlabel='$\Delta R$(Z/H, 2nd closest b) (Zh_pt>300)', bin_range=[0,3], bins=15, add_cuts='Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_2ndcloseb_dr',xlabel='$\Delta R$(Z/H, 2nd closest b) (Zh_pt>450)', bin_range=[0,3], bins=15, add_cuts='Zh_pt>450', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_closeq_dr',xlabel='$\Delta R$(Z/H, closest q) ',          bin_range=[0,3], bins=15,  sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_closeq_dr',xlabel='$\Delta R$(Z/H, closest q) (Zh_pt>300)', bin_range=[0,3], bins=15, add_cuts='Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        #Hist(['ttZ','ttH'],'Zh_closeq_dr',xlabel='$\Delta R$(Z/H, closest q) (Zh_pt>450)', bin_range=[0,3], bins=15, add_cuts='Zh_pt>450', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
        ## --- End of study for andrew
        
        ## --- Study for BTV comments
        #alt_weight_1 = (lambda v_, obj: v_['weight']* np.sign(v_['genWeight'])* v_['topptWeight']
        #                * (v_['HEM_weight'] if obj.year+obj.HEM_opt == '2018' else 1.0)* (v_['lep_trigeffsf'])* v_['lep_sf']
        #                * v_['dak8md_bbvl_sf']* v_['puWeight']* (v_['PrefireWeight'] if obj.year != '2018' else 1.0)
        #)
        #alt_weight_2 = (lambda v_, obj: v_['weight']* np.sign(v_['genWeight'])* v_['topptWeight']
        #                * (v_['HEM_weight'] if obj.year+obj.HEM_opt == '2018' else 1.0)* (v_['lep_trigeffsf'])* v_['lep_sf']
        #                * v_['dak8md_bbvl_sf']* v_['puWeight']* (v_['PrefireWeight'] if obj.year != '2018' else 1.0)
        #                * v_['BTagWeight_nocorr']
        #)
        #alt_weight_3 = (lambda v_, obj: v_['weight']* np.sign(v_['genWeight'])* v_['topptWeight']
        #                * (v_['HEM_weight'] if obj.year+obj.HEM_opt == '2018' else 1.0)* (v_['lep_trigeffsf'])* v_['lep_sf']
        #                * v_['dak8md_bbvl_sf']* v_['puWeight']* (v_['PrefireWeight'] if obj.year != '2018' else 1.0)
        #                * v_['BTagWeight']
        #)
        #
        #StackedHist(processes,'HT', xlabel='HT (no btagSF)',         bin_range=[0,2000],  n_bins=40,  doCuts=True, addData=True, doShow=False, alt_weight=alt_weight_1)  
        #StackedHist(processes,'HT', xlabel='HT (with btagSF+corr)',  bin_range=[0,2000],  n_bins=40,  doCuts=True, addData=True, doShow=False, alt_weight=alt_weight_3)  
        ##
        #StackedHist(processes, 'jetbtag_1', xlabel=r'leading jet deepCSV (AK4, with btagSF)',       bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False, alt_weight=alt_weight_2)  
        #StackedHist(processes, 'jetbtag_1', xlabel=r'leading jet deepCSV (AK4, with btagSF+corr) ', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False, alt_weight=alt_weight_3)  
        ##
        #StackedHist(processes, 'n_ak4jets', xlabel='# of ak4 jets (no btagSF)',      bin_range=[-0.5,12.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5],  doCuts=True, addData=True, doShow=False, alt_weight=alt_weight_1) 
        #StackedHist(processes, 'n_ak4jets', xlabel='# of ak4 jets (with btagSF+corr)',  bin_range=[-0.5,12.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5],  doCuts=True, addData=True, doShow=False, alt_weight=alt_weight_3) 


        ## --- end of BTV study

        ## --- BBvL Fakse SF control region study
        #Plotter.reset_cut_func()
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'n_genb_matchZH', xlabel='n_genb_matchZH (Analysis SR, ZH bbvl > 0.8) ', bins=[-0.5,0.5,1.5,2.5,3.5],  add_cuts='Zh_bbvLscore>0.8',   doCuts=True, doLog=False, doNorm=True, doShow=False)
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'NN', xlabel='NN (Analysis SR, ZH bbvl > 0.8, n_genb_matched = 0) ', bin_range=[0,1], n_bins=20,  add_cuts='Zh_bbvLscore>0.8;n_genb_matchZH==0',   doCuts=True, doLog=True, doNorm=True, doShow=False)
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'NN', xlabel='NN (Analysis SR, ZH bbvl > 0.8, n_genb_matched > 0) ', bin_range=[0,1], n_bins=20,  add_cuts='Zh_bbvLscore>0.8;n_genb_matchZH>=1',   doCuts=True, doLog=True, doNorm=True, doShow=False)
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'Zh_M', xlabel='Z/H M (Analysis SR, ZH bbvl > 0.8, n_genb_matched = 0) ', bin_range=[50,200], n_bins=15,  add_cuts='Zh_bbvLscore>0.8;n_genb_matchZH==0',   doCuts=True, doLog=False, doNorm=True, doShow=False)
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'Zh_M', xlabel='Z/H M (Analysis SR, ZH bbvl > 0.8, n_genb_matched > 0) ', bin_range=[50,200], n_bins=15,  add_cuts='Zh_bbvLscore>0.8;n_genb_matchZH>=1',   doCuts=True, doLog=False, doNorm=True, doShow=False)
        #Plotter.set_cut_func(getFakebbvlCuts)
        #alt_weight = getFakebbvlWeights
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'n_genb_matchZH', xlabel='n_genb_matchZH (Analysis CR, ZH bbvl > 0.8) ', bins=[-0.5,0.5,1.5,2.5,3.5],  add_cuts='Zh_bbvLscore>0.8',   alt_weight=alt_weight, doCuts=True, doLog=False, doNorm=True, doShow=False)
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'Zh_bbvLscore', xlabel='deepAK8MD bbvL (Analysis CR, n_genb_matched = 0) ', bin_range=[0,1], n_bins=20,  add_cuts='n_genb_matchZH==0',   alt_weight=alt_weight, doCuts=True, doLog=False, doNorm=True, doShow=False)
        #Hist(['TTBar','tt_B','ttZ','ttH'], 'Zh_bbvLscore', xlabel='deepAK8MD bbvL (Analysis CR, n_genb_matched > 0) ', bin_range=[0,1], n_bins=20,  add_cuts='n_genb_matchZH>=1',   alt_weight=alt_weight, doCuts=True, doLog=False, doNorm=True, doShow=False)

        #
        #StackedHist(processes,    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV]', bin_range=[200,600], n_bins=20, alt_weight=alt_weight,  doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_eta', xlabel=r'Z/H $\eta$', bin_range=[-2.5,2.5], n_bins=20, alt_weight=alt_weight, doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$', bin_range=[-3.15,3.15], n_bins=20, alt_weight=alt_weight, doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,   'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV]', bin_range=[50,200],  n_bins=15, alt_weight=alt_weight, doCuts=True,  doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bbvLscore', xlabel='Z/H deepAK8MD bbvL', bin_range=[0,1],  n_bins=20, alt_weight=alt_weight,  doCuts=True, addData=True, doLog=True, doShow=False)  
        #
        ##
        #StackedHist(processes,    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV] (bbvL > 0.8)', bin_range=[200,600], n_bins=20, alt_weight=alt_weight, add_cuts='Zh_bbvLscore>0.8',  doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_eta', xlabel=r'Z/H $\eta$ (bbvL > 0.8)', bin_range=[-2.5,2.5], n_bins=20,  doCuts=True, alt_weight=alt_weight, add_cuts='Zh_bbvLscore>0.8', doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$ (bbvL > 0.8)', bin_range=[-3.15,3.15], n_bins=20,  doCuts=True, alt_weight=alt_weight, add_cuts='Zh_bbvLscore>0.8', doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,   'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV] (bbvL > 0.8)', bin_range=[50,200],  n_bins=15, alt_weight=alt_weight, add_cuts='Zh_bbvLscore>0.8', doCuts=True,  doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bbvLscore', xlabel='Z/H deepAK8MD bbvL', bins=[0.8,.97,1], add_cuts='Zh_bbvLscore>0.8', alt_weight=alt_weight, doCuts=True, addData=True, doLog=True, doShow=False)  

        ## --- end of bbvl fake sf CR study


        #
        #Hist(['ttZ','ttH'],'withbbvl_NN',xlabel='NN, matchedGen_ZHbb;Zh_pt>450', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80;Zh_pt>450', doNorm=False, doCuts=True, doShow=False)  
        ###StackedHist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='NN>0.80;Zh_pt>300', doCuts=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withdak8md_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withdak8md_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (noak8md_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='noak8md_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets'],'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withbbvl_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doBinRatio=True, doLog=False, addData=False, doShow=False)  
        #Hist(['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets'],'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,75,105,145,200], add_cuts='withbbvl_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doBinRatio=True, doLog=False, addData=False, doShow=False)  
        #Hist(['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets'],'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,150,200], add_cuts='withbbvl_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doBinRatio=True, doLog=False, addData=False, doShow=False)  
        #Hist(['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets'],'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,75,105,150,200], add_cuts='withbbvl_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doBinRatio=True, doLog=False, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (ZhpT>300)', bin_range=[50,200],  bins=30, add_cuts='Zh_pt>300', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  

        # ---- ana strat

        #StackedHist(processes,    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV]', bins=[200,300,450,600],  doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,'NN', xlabel='DNN score', bin_range=[0,1],  n_bins=20,  add_d_cuts='NN<=0.8',     doCuts=True, addData=True, doShow=False)  
        #Hist(['ttZ','ttH','TTBar','tt_B'],'NN', xlabel='DNN score', bin_range=[0,1],  n_bins=10, doNorm=True, doLog=True, doCuts=True, addData=False, doShow=False)  
        #StackedHist(processes,'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV]', bin_range=[50,200],  bins=cfg.sdm_bins, doCuts=True,  doLog=True, addData=True, doShow=False)  
        #
        #Hist(['ttZ','ttH','TTBar','tt_B'],'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV] (Z/H $p_{\mathrm{T}}>300$ (GeV);DNN$>0.8$)', bin_range=[50,200],  bins=[50,75,90,105,120,140,200], add_cuts='Zh_pt>300;NN>0.80', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  
        #Hist(['ttZ','ttH','TTBar','tt_B'],'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV] (Z/H $p_{\mathrm{T}}>450$ (GeV);DNN$>0.8$)', bin_range=[50,200],  bins=[50,75,90,105,120,140,200], add_cuts='Zh_pt>450;NN>0.80', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  

        # ---- end ana strat

        # ---- Zh info

        #StackedHist(processes,    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV]', bin_range=[200,600], n_bins=20,  doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_eta', xlabel=r'Z/H $\eta$', bin_range=[-2.5,2.5], n_bins=20,  doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$', bin_range=[-3.15,3.15], n_bins=20,  doCuts=True, doLog=True, addData=True, doShow=False)  
        #StackedHist(processes,   'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV]', bin_range=[50,200],  n_bins=15, doCuts=True,  doLog=True, addData=True, doShow=False)  
        #
        #StackedHist(processes,    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ (Zhpt>300, NN>0.35) [GeV]', bin_range=[200,600], n_bins=20, add_cuts='Zh_pt>300;NN>0.35',  doCuts=True, doLog=False, addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_eta', xlabel=r'Z/H $\eta$ (Zhpt>300, NN>0.35)', bin_range=[-2.5,2.5], n_bins=20,  doCuts=True, doLog=False, add_cuts='Zh_pt>300;NN>0.35', addData=True, doShow=False)  
        #StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$ (Zhpt>300, NN>0.35)', bin_range=[-3.15,3.15], n_bins=20,  doCuts=True, doLog=False, add_cuts='Zh_pt>300;NN>0.35', addData=True, doShow=False)  
        #StackedHist(processes,   'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ (Zhpt>300, NN>0.35) [GeV]', bin_range=[50,200],  n_bins=15, doCuts=True, add_cuts='Zh_pt>300;NN>0.35',  doLog=False, addData=True, doShow=False)  

        # ---- end zh info

        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (ZhpT>450)', bin_range=[50,200],  bins=30, add_cuts='Zh_pt>450', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (ZhpT>450;NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='Zh_pt>450;withbbvl_NN>0.80', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  
        ##StackedHist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (NN>0.90,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='NN>0.90;Zh_pt>300', doCuts=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withdak8md_NN>0.90,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withdak8md_NN>0.90;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (noak8md_NN>0.90,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='noak8md_NN>0.90;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.90,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withbbvl_NN>0.90;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        ##StackedHist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (NN>0.95,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='NN>0.95;Zh_pt>300', doCuts=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withdak8md_NN>0.95,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withdak8md_NN>0.95;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (noak8md_NN>0.95,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='noak8md_NN>0.95;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.95,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withbbvl_NN>0.95;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  

        #StackedHist(processes,'Zh_doubleB', xlabel='Zh_doubleB (NN>.90, 80<Zh_M<145)', bin_range=[-1,1],  n_bins=20, add_cuts='NN>=.90;Zh_M>80;Zh_M<145', doCuts=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_doubleB', xlabel='Zh_doubleB (NN>.90, 80<Zh_M<145)', bin_range=[-1,1],  n_bins=20, add_cuts='NN>=.90;Zh_M>80;Zh_M<145', doCuts=True, doNorm=False, doLog=True, doShow=False)  
        #Hist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLsscore (n_b_outZH>=2)', bin_range=[0,1], add_cuts='n_b_outZh>=2',  n_bins=20, doCuts=False, doNorm=True, doShow=False) 
        # ----
        #StackedHist(processes,'b1_outZh_score', bin_range=[.4,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_deepB', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bbvLscore', bin_range=[0,1],  n_bins=20,  doCuts=True, add_cuts='Zh_bbvLscore>0.6', addData=True, doShow=False, sepGenOpt='sepGenSig')  
        #StackedHist(processes,'best_rt_score', bin_range=[0,1],  n_bins=20,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'ak4_worstb_inZH', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=True)  
        #StackedHist(processes,'ak4_bestb_inZH',  bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=True)  
        #StackedHist(processes,'Zh_bbscore_sj', bin_range=[0,2],  n_bins=20,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_Tscore', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_Wscore', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        
        # --- sdm for ken
        
        #Hist(['ttH','ttZ'],    'fjetsdm_1', xlabel=r'leading fatjet $m_{sd}$ [GeV] (gen matched)', bin_range=[50,200], add_cuts='matchedGen_ZHbb_bb==1',  n_bins=30,  doCuts=False, doNorm=False, doShow=False)  

        # --- end sdm
        #StackedHist(['ttH','ttZ','tt_B','TTBar'],    'NN',  n_bins=20, bin_range=[0,1], doNorm=False,  doCuts=True, doLog=True, doShow=True)      
        #Hist(['ttH','ttZ'],    'Zh_pt', xlabel=r'RECO Z/H $p_{T}$ [GeV] (450 cut)', n_bins=20, bin_range=[0,1000], doNorm=False, add_cuts='Zh_pt>450',  doCuts=True, doLog=True, doShow=False)      
        #Hist(['ttH','ttZ'],    'genZHpt', xlabel=r'GEN Z/H $p_{T}$ [GeV]', n_bins=20, bin_range=[0,1000], doNorm=False,  doCuts=True, doLog=True, doShow=False)      
        #Hist(['ttH','ttZ'],    'genZHpt', xlabel=r'GEN Z/H $p_{T}$ [GeV]', bins=[0,100,200,300,450,700,1200], doNorm=False,  doCuts=True, doLog=True, doShow=False)  


        # --- trigger
        
        #StackedHist(processes,    'pbt_muon', xlabel= r'Muon trigger', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        #StackedHist(processes,    'HLT_Mu50', xlabel= r'HLT_Mu50', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        #if y != '2017':
        #    StackedHist(processes,    'HLT_IsoMu24', xlabel= r'HLT_IsoMu24', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        #else:
        #    StackedHist(processes,    'HLT_IsoMu27', xlabel= r'HLT_IsoMu27', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        ##
        #if y == '2016':
        #    StackedHist(processes,    'HLT_IsoTkMu24', xlabel= r'HLT_IsoTkMu24', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        #    StackedHist(processes,    'HLT_TkMu50', xlabel= r'HLT_TkMu50', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        #if int(y) >= 2017: 
        #    StackedHist(processes,    'HLT_OldMu100', xlabel= r'HLT_OldMu100', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  
        #    StackedHist(processes,    'HLT_TkMu100', xlabel= r'HLT_TkMu100', bins=[-.5,0.5,1.5],   add_cuts='passSingleLepMu==1;MET_pt>20;pbt_muon==1', doCuts=False, addData=True,  doShow=False)  

        Hist(['tt_B'],    'LHE_HTIncoming', xlabel='LHE_HTIncoming (inc)', n_bins=50, bin_range=[0, 1500], doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
        Hist(['tt_B'],    'LHE_HTIncoming', xlabel='LHE_HTIncoming (300 < ZHpt < 450)', n_bins=50, bin_range=[0, 1500], add_cuts='Zh_pt>300;Zh_pt<450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
        Hist(['tt_B'],    'LHE_HTIncoming', xlabel='LHE_HTIncoming (ZHpt > 450)', n_bins=50, bin_range=[0, 1500],       add_cuts='Zh_pt>450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  

        Hist(['tt_B'],    'LHE_HT', xlabel='LHE_HT (inc)', n_bins=50, bin_range=[0, 1500], doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
        Hist(['tt_B'],    'LHE_HT', xlabel='LHE_HT (300 < ZHpt < 450)', n_bins=50, bin_range=[0, 1500], add_cuts='Zh_pt>300;Zh_pt<450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
        Hist(['tt_B'],    'LHE_HT', xlabel='LHE_HT (ZHpt > 450)', n_bins=50, bin_range=[0, 1500],       add_cuts='Zh_pt>450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  


    #
    return 1



if __name__ == '__main__':
    main()
