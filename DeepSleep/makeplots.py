from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
from modules.AnaDict import AnaDict
import cfg.deepsleepcfg as cfg
###----=======-----###

#sepGenOpt ; 'sepGenSig','sepGenBkg','sepGenMatchedSig','sepGenMatchedBkg'

###----=======-----###


#Plotter.load_data('2018','postHEM')
#Plotter.load_data('2017')
#StackedHist(cfg.MC_samples,'Lep_pt', xlabel='Elec_pt',   bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)


for y in cfg.Years:
    Plotter.load_data(y)

    #Hist(['TTZH','TTBarLep'],'NN', bin_range=[0,1], doNorm=True,  n_bins=20, sepGenOpt='sepGenSig;sepGenBkg', doCuts=True, addData=False)
    StackedHist(cfg.MC_samples,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=0.8', doCuts=True, addData=True)

    eff_tag  = 'tight'

    mu_eff     = AnaDict.read_pickle(f'files/{y}/eff_files/muon_eff_{y}_{eff_tag}.pkl')
    ele_eff    = AnaDict.read_pickle(f'files/{y}/eff_files/electron_eff_{y}_{eff_tag}.pkl')
    #
    mu_pt_bins = [float(bin_str.split(',')[0]) for bin_str in mu_eff['pt']] + [float(list(mu_eff['pt'].keys())[-1].split(',')[1])]
    ele_pt_bins = [float(bin_str.split(',')[0]) for bin_str in ele_eff['pt']] + [float(list(ele_eff['pt'].keys())[-1].split(',')[1])]
    mu_eta_bins = [float(bin_str.split(',')[0]) for bin_str in mu_eff['eta']] + [float(list(mu_eff['eta'].keys())[-1].split(',')[1])]
    ele_eta_bins = [float(bin_str.split(',')[0]) for bin_str in ele_eff['eta']] + [float(list(ele_eff['eta'].keys())[-1].split(',')[1])]

    
    
    
    #StackedHist(cfg.MC_samples,'Lep_pt', xlabel='Elec_pt_pass_trigger_electron_pbt (GeV)',   bin_range=[0,700],      bins=ele_pt_bins,    add_cuts='MET_pt>=20;passSingleLepElec==1;pbt_elec==1',    doCuts=False, addData=True)
    #StackedHist(cfg.MC_samples,'Lep_eta', xlabel='Elec_eta_pass_trigger_electron_pbt', bin_range=[-2.6,2.6],   bins=ele_eta_bins,   add_cuts='MET_pt>=20;passSingleLepElec==1;pbt_elec==1',    doCuts=False, addData=True)
    #
    #
    #StackedHist(cfg.MC_samples,'Lep_pt',   xlabel='Mu_pt_pass_trigger_muon_pbt (GeV)',  bin_range=[0,700],    bins=mu_pt_bins,   add_cuts='MET_pt>=20;passSingleLepMu==1;pbt_muon==1',  doCuts=False, addData=True)
    #StackedHist(cfg.MC_samples,'Lep_eta',  xlabel='Mu_eta_pass_trigger_muon_pbt', bin_range=[-2.6,2.6],  bins=mu_eta_bins,   add_cuts='MET_pt>=20;passSingleLepMu==1;pbt_muon==1',  doCuts=False, addData=True)

    
exit()
#StackedHist(cfg.MC_samples,'Pass_trigger_electron',    bin_range=[-1,2],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=3, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Pass_trigger_muon',    bin_range=[-1,2],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=3, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=18, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phiM',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_eta<0.0',    n_bins=18, doCuts=False, addData=True)

#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi_ptg200',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_pt>200',    n_bins=18, doCuts=False, addData=True)

#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt preHEM',  bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt preHEM',   bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_eta',  xlabel='Ele_eta_preHEM',   bin_range=[-3,3],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi_EC_ptg200_preHEM',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_pt>200;Lep_eta<-1.4',    n_bins=18, doCuts=False, addData=True)
Plotter.load_data('2018','postHEM')                 
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt postHEM',  bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1;SAT_Pass_HEMVeto_DataAndMC==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt postHEM',   bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1;SAT_Pass_HEMVeto_DataAndMC==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_eta',  xlabel='Ele_eta_postHEM',   bin_range=[-3,3],   add_cuts='MET_pt>=20;passSingleLepElec==1;SAT_Pass_HEMVeto_DataAndMC==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi_EC_ptg200_postHEM',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_pt>200;Lep_eta<-1.4;SAT_Pass_HEMVeto_DataAndMC==1',    n_bins=18, doCuts=False, addData=True)
#
#StackedHist(cfg.MC_samples,'Lep_pt',    bin_range=[0,500],   add_cuts='MET_pt>=20',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'MET_pt',    bin_range=[0,500],   add_cuts='MET_pt>=20',   n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'n_ak4jets', bin_range=[2.5,14.5],      add_cuts='MET_pt>=20',   n_bins=12, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'n_ak8jets', bin_range=[-0.5,6.5],       add_cuts='MET_pt>=20',   n_bins=7, doCuts=False, addData=True)
#
#
#
#Plotter.load_data('2016')
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt',  bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt',   bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',    bin_range=[0.0,500],   add_cuts='MET_pt>=20',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'MET_pt',    bin_range=[0.0,500],   add_cuts='MET_pt>=20',   n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'n_ak4jets', bin_range=[2,14],      add_cuts='MET_pt>=20',   n_bins=14, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'n_ak8jets', bin_range=[0,6],       add_cuts='MET_pt>=20',   n_bins=8, doCuts=False, addData=True)

Plotter.load_data('2017')
#StackedHist(cfg.MC_samples,'Pass_trigger_electron',     bin_range=[-1,3],   add_cuts='MET_pt>=20;passSingleLepElec==1', n_bins=10, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt',  bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt',   bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',    bin_range=[0.0,500],   add_cuts='MET_pt>=20',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'MET_pt',    bin_range=[0.0,500],   add_cuts='MET_pt>=20',   n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'n_ak4jets', bin_range=[2.5,14.5],      add_cuts='MET_pt>=20',   n_bins=12, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'n_ak8jets', bin_range=[-0.5,6.5],       add_cuts='MET_pt>=20',   n_bins=7, doCuts=False, addData=True)

