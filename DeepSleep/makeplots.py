from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
import cfg.deepsleepcfg as cfg
from functools import reduce
###----=======-----###

#sepGenOpt :'sepGenSig','sepGenBkg','sepGenMatchedSig','sepGenMatchedBkg'

###----=======-----###

#StackedHist(cfg.MC_samples,'Pass_trigger_electron',     bin_range=[-1,3],   add_cuts='MET_pt>=20', n_bins=10, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Pass_trigger_muon',     bin_range=[-1,3],   add_cuts='MET_pt>=20', n_bins=10, doCuts=False, addData=True)

Plotter.load_data('2018','preHEM')
#StackedHist(cfg.MC_samples,'SAT_Pass_HEMVeto_DataOnly', xlabel='SAT_Pass_HEMVeto_Data_mu', bin_range=[-1,2],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=4, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'SAT_Pass_HEMVeto_DataOnly', xlabel='SAT_Pass_HEMVeto_DataA_ele', bin_range=[-1,2],   add_cuts='MET_pt>=20;passSingleLepElec==1;Pass_trigger_electron==1',    n_bins=4, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_eta',  xlabel='Ele_eta',   bin_range=[-3,3],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=18, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phiM',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_eta<0.0',    n_bins=18, doCuts=False, addData=True)

#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi_ptg200',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_pt>200',    n_bins=18, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_phi',  xlabel='Ele_phi_EC_ptg200',   bin_range=[-3.14,3.14],   add_cuts='MET_pt>=20;passSingleLepElec==1;Lep_pt>200;Lep_eta<-1.4',    n_bins=18, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt preHEM',  bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt preHEM',   bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
Plotter.load_data('2018','postHEM')                 
StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt postHEM',  bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1;SAT_Pass_HEMVeto_DataAndMC==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt postHEM',   bin_range=[0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1;SAT_Pass_HEMVeto_DataAndMC==1',    n_bins=20, doCuts=False, addData=True)
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

