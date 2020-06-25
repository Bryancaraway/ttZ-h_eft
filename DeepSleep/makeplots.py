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
#Plotter.load_data('2016')
#StackedHist(cfg.MC_samples,'Lep_pt',    xlabel='Ele pt',  bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1', n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Lep_pt',    xlabel='Mu pt',   bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1',   n_bins=20, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'MET_pt', bin_range=[0.0,500],   add_cuts='MET_pt>=20',   n_bins=20, doCuts=False, addData=True)

Plotter.load_data('2017')
StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Ele_pt',  bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepMu==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_pt',  xlabel='Mu_pt',   bin_range=[0.0,500],   add_cuts='MET_pt>=20;passSingleLepElec==1',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_pt',    bin_range=[0.0,500],   add_cuts='MET_pt>=20',    n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'MET_pt',    bin_range=[0.0,500],   add_cuts='MET_pt>=20',   n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'n_ak4jets', bin_range=[2,14],      add_cuts='MET_pt>=20',   n_bins=14, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'n_ak8jets', bin_range=[0,6],       add_cuts='MET_pt>=20',   n_bins=8, doCuts=False, addData=True)

