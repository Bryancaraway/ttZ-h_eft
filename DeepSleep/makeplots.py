from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
import cfg.deepsleepcfg as cfg
from functools import reduce


#sepGenOpt :'sepGenSig','sepGenBkg','sepGenMatchedSig','sepGenMatchedBkg'

Plotter.load_data('2016')


#StackedHist(cfg.MC_samples,'NN', bin_range=[0.0,1.0],
#            n_bins=20,
#            #sepGenOpt='sepGenSig;sepGenBkg',
#            add_cuts='NN>=0.0;NN<=0.8',
#            addData=True)
#StackedHist(cfg.MC_samples,'Pass_trigger_electron',     bin_range=[-1,3],   add_cuts='MET_pt>=20', n_bins=10, doCuts=False, addData=True)
#StackedHist(cfg.MC_samples,'Pass_trigger_muon',     bin_range=[-1,3],   add_cuts='MET_pt>=20', n_bins=10, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'weight',    xlabel='weights',  bin_range=[-5,5], n_bins=100, doCuts=False, addData=False)
Plotter.load_data('2017')
StackedHist(cfg.MC_samples,'weight',    xlabel='weights',  bin_range=[-5,5], n_bins=100, doCuts=False, addData=False)
StackedHist(cfg.MC_samples,'Lep_pt',    xlabel='Ele pt',  bin_range=[0.0,500],   add_cuts='MET_pt>=20;Pass_trigger_electron==1', n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'Lep_pt',    xlabel='Mu pt',   bin_range=[0.0,500],   add_cuts='MET_pt>=20;Pass_trigger_muon==1',     n_bins=20, doCuts=False, addData=True)

StackedHist(cfg.MC_samples,'Lep_eta',     bin_range=[-2.6,2.6],  add_cuts='MET_pt>=20', n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'MET_pt',      bin_range=[0,500],     add_cuts='MET_pt>=20', n_bins=20, doCuts=False, addData=True)
StackedHist(cfg.MC_samples,'spher',       bin_range=[0,1],       add_cuts='MET_pt>=20', n_bins=10, doCuts=True,  addData=True)
StackedHist(cfg.MC_samples,'n_ak4jets',   bin_range=[0,12],      add_cuts='MET_pt>=20', n_bins=12, doCuts=False,  addData=True)
StackedHist(cfg.MC_samples,'n_ak8jets',   bin_range=[0,5],       add_cuts='MET_pt>=20', n_bins=5,  doCuts=False,  addData=True)
