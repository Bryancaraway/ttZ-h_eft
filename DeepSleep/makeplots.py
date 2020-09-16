
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
from modules.AnaDict import AnaDict
import config.ana_cff as cfg
###----=======-----###

#sepGenOpt ; 'sepGenSig','sepGenBkg','sepGenMatchedSig','sepGenMatchedBkg'

###----=======-----###

#jjec = 'ak8'
#jec_list = ['JESUp','JESDown','JERUp','JERDown']


for y in cfg.Years: 
    #for jec in jec_list:
    #Plotter.load_data(y, addBSF=False, tag=f'{jjec}{jec}') #tag='ak4JESUp'
    Plotter.load_data(y, addBSF=False)
    ''' LOOK AT STACKED DATA VS MC '''
    #StackedHist(cfg.MC_pow,    'Zh_pt', bin_range=[200,500],  n_bins=20,  doCuts=False, addData=True)  
        #StackedHist(cfg.MC_pow,    'Zh_pt', bin_range=[200,500],  n_bins=20,  doCuts=True, addData=True)  
        #
        #StackedHist(cfg.MC_pow,    'Zh_M', bin_range=[0,250],  bins=[50,80,105,145,200],  doCuts=False, addData=True)  
        #StackedHist(cfg.MC_pow,    'Zh_M', bin_range=[0,250],  n_bins=20,  doCuts=True, addData=True)  
        #
        #StackedHist(cfg.MC_pow,    'n_b_outZh', bin_range=[-.5,3.5],  bins=[-.5,.5,1.5,2.5,3.5],  doCuts=False, addData=True)  
        #StackedHist(cfg.MC_pow,    'n_b_outZh', bin_range=[-.5,3.5],  bins=[-.5,.5,1.5,2.5,3.5],  doCuts=False, add_cuts='n_b_outZh>1', addData=True)  
        #StackedHist(cfg.MC_pow,    'n_b_outZh', bin_range=[-.5,3.5],  bins=[-.5,.5,1.5,2.5,3.5],  doCuts=True, addData=True)  
        #
        #StackedHist(cfg.MC_pow,    'n_q_outZh', bin_range=[-.5,3.5],  bins=[-.5,.5,1.5,2.5,3.5],  doCuts=False, addData=True)  
        #StackedHist(cfg.MC_pow,    'n_q_outZh', bin_range=[-.5,3.5],  bins=[-.5,.5,1.5,2.5,3.5],  doCuts=True, addData=True)  
        ##
        #StackedHist(cfg.MC_pow,    'nonZhbb_b1_dr', bin_range=[0,5.5],  n_bins=10,  doCuts=False, addData=True)  
        #StackedHist(cfg.MC_pow,    'nonZhbb_b1_dr', bin_range=[0,5.5],  n_bins=10,  doCuts=True, addData=True)  
        ##
        #StackedHist(cfg.MC_pow,    'nonZhbb_b2_dr', bin_range=[0,5.5],  n_bins=10,  doCuts=False, addData=True)  
        #StackedHist(cfg.MC_pow,    'nonZhbb_b2_dr', bin_range=[0,5.5],  n_bins=10,  doCuts=True, addData=True)  
        ####====####

    #StackedHist(cfg.MC_pow,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='NN<=0.8', doCuts=True, addData=True)  
    StackedHist(cfg.MC_pow,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='NN<=0.8', sepGenOpt='sepGenBkg;--', doCuts=True, addData=True)  

    #StackedHist(cfg.MC_pow,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=0.80', doCuts=True, sepGenOpt='sepGenBkg;++', addData=True)  
    #StackedHist(cfg.MC_pow,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', doCuts=True, sepGenOpt='sepGenBkg;++', addData=False)  
    #StackedHist(cfg.MC_pow,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=0.80', doCuts=True, addData=True)  
    #StackedHist(cfg.MC_samples,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=.80', doCuts=True, addData=True)  
    #StackedHist(cfg.MC_pow,    'n_ak4jets', bin_range=[4,12],  bins=[4,5,6,7,8,9,10,11,12],  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'n_ak4jets', bin_range=[4,12],  bins=[4,5,6,7,8,9,10,11,12],  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'n_ak8jets', bin_range=[0,4],  bins=[0,1,2,3,4],  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'n_ak8jets', bin_range=[0,4],  bins=[0,1,2,3,4],  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'nBottoms_drLeptonCleaned', bin_range=[-0.5,8.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_pow,    'nBottoms_drLeptonCleaned', bin_range=[-0.5,8.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],  doCuts=True, addData=True)  
    #
    #StackedHist(cfg.MC_pow,    'nBottoms_drLeptonCleaned', bin_range=[0,8],  bins=[0,1,2,3,4,5,6,7,8],  doCuts=True, addData=True)  
    #StackedHist(cfg.MC_samples,'nBottoms_drLeptonCleaned', bin_range=[0,8],  bins=[0,1,2,3,4,5,6,7,8],  doCuts=False, addData=True)  
    #
    #StackedHist(cfg.MC_pow,    'spher', bin_range=[0,1],  n_bins=10,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'spher', bin_range=[0,1],  n_bins=10,  doCuts=False, addData=True)  
    #
    #StackedHist(cfg.MC_pow,    'MET_pt', bin_range=[0,300],  n_bins=15,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'MET_pt', bin_range=[0,300],  n_bins=15,  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'Zh_M', bin_range=[50,200],  n_bins=15,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_pow,    'Zh_M', bin_range=[50,200],  n_bins=15,  doCuts=True, addData=True)  
    #StackedHist(cfg.MC_samples,'Zh_M', bin_range=[50,200],  n_bins=15,  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'Zh_pt', bin_range=[200,500],  n_bins=20,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'Zh_pt', bin_range=[200,500],  n_bins=20,  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'Zh_score', bin_range=[0,1],  n_bins=10,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'Zh_score', bin_range=[-1,1],  n_bins=20,  doCuts=False, addData=True)  
    #
    #StackedHist(cfg.MC_pow,    'n_b_inZh', bin_range=[-0.5,2.5],  n_bins=3,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'n_b_inZh', bin_range=[-0.5,2.5],  n_bins=3,  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'n_b_outZh', bin_range=[-0.5,4.5],  n_bins=5,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'n_b_outZh', bin_range=[-0.5,4.5],  n_bins=5,  doCuts=False, addData=True)  
    ##
    #StackedHist(cfg.MC_pow,    'n_ak8_Zhbb', bin_range=[-0.5,4.5],  n_bins=5,  doCuts=False, addData=True)  
    #StackedHist(cfg.MC_samples,'n_ak8_Zhbb', bin_range=[-0.5,4.5],  n_bins=5,  doCuts=False, addData=True)  
    
    
    ''' STUDY MG vs POWHEG '''
    
    
    #Hist(['TTBarLep','TTBarLep_pow'],'NN', bin_range=[0,1],       doNorm=False, doLog=True, n_bins=20,  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow'],'Zh_pt', bin_range=[200,500],  doNorm=False, doLog=True, n_bins=20,  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow'],'Zh_M', xlabel='Zh_M (NN>0.9;Zh_pt>300)', bin_range=[50,200], bins=[50,80,105,145,200],  doNorm=False, doLog=False, add_cuts='NN>0.9;Zh_pt>300',  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow'],'nBottoms_drLeptonCleaned', bin_range=[0,7] ,bins=[0,1,2,3,4,5,6,7],  doNorm=False, doLog=True,  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow'],'Zh_score', bin_range=[0,1], n_bins=20,  doNorm=False, doLog=True,  doCuts=True, addData=False)
    
    ''' STUDY MG tt+bb VS POWHEG tt+bb VS POWHEG DEDICATED tt+bb '''
    
    #Hist(['TTBarLep','TTBarLep_pow','TT_bb_pow'],'Zh_M', xlabel='Zh_M (NN>0.9;Zh_pt>300)', bin_range=[50,200], bins=[50,80,105,145,200], doNorm=True,  sepGenOpt='sepGenBkg', droptt=True,   add_cuts='NN>0.9;Zh_pt>300',  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow','TT_bb_pow'],'Zh_pt', bin_range=[200,500], doNorm=True,  sepGenOpt='sepGenBkg', droptt=True, n_bins=20, doCuts=True, addData=False)
    #Hist(['TTBarLep_pow','TTZ_bb','TTZH','TT_bb_pow'],'Zh_M', xlabel='Zh_M (NN>0.9;Zh_pt>300)', bin_range=[50,200],bins=[50,80,105,145,200], doNorm=True,  sepGenOpt='sepGenSig;sepGenBkg', dropZqq=True,  add_cuts='NN>0.9;Zh_pt>300',  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow','TT_bb_pow'],'Zh_score', bin_range=[0,1], doNorm=True,  sepGenOpt='sepGenBkg', droptt=True, n_bins=20,  doCuts=True, addData=False)
    #Hist(['TTBarLep','TTBarLep_pow','TT_bb_pow'],'NN', xlabel='NN_tt+bb_pow', bin_range=[0,1], doNorm=True,  sepGenOpt='sepGenBkg', droptt=True, n_bins=20,  doCuts=True, addData=False)
    
    ''' STUDY TTZ,Zbb vs DEDICATED TTZ,Zbb '''
    
    #Hist(['TTZH','TTBarLep_pow'], 'Zh_M', xlabel='Zh_M (NN>0.9;Zh_pt>300)', bin_range=[50,200], bins=[50,80,105,145,200], doNorm=True,  sepGenOpt='sepGenSig;sepGenBkg;++', droptt=True, dropZqq=True, add_cuts='NN>0.9;Zh_pt>300', doCuts=True, addData=False)
    #Hist(['TTZ_bb','TTZH'],'Zh_M', xlabel='Zh_M (NN>0.9;Zh_pt>300)', bin_range=[50,200],bins=[50,80,105,145,200], doNorm=True,  sepGenOpt='sepGenSig', dropZqq=True,  add_cuts='NN>0.9;Zh_pt>300', doCuts=True, addData=False)
    #Hist(['TTZH'],'Zh_M', xlabel='Zh_M (NN>0.9;Zh_pt>300)', bin_range=[50,200],bins=[50,80,105,145,200], doNorm=True,  sepGenOpt='sepGenSig', dropZqq=True,  add_cuts='NN>0.9;Zh_pt>300', doCuts=True, addData=False)
    #Hist(['TTZH'],'NN', bin_range=[0,1], doNorm=True,  sepGenOpt='sepGenSig',  n_bins=20, dropZqq=True,  doCuts=True, addData=False)
    
    #Hist(['TTZH','TTZ_bb'],'n_b_outZh', bin_range=[-0.5,4.5], n_bins=5, doNorm=True,  sepGenOpt='sepGenSig', dropZqq=True,  doCuts=False, addData=False)
    #Hist(['TTZH','TTZ_bb'],'Zh_score', bin_range=[0,1], doNorm=True,  sepGenOpt='sepGenSig',  n_bins=20, dropZqq=True,  doCuts=True, addData=False)
exit()
#

for y in ['2018']:#cfg.Years:
    Plotter.load_data(y)
    

    #Hist(['TTZH','TTBarLep'],'NN', bin_range=[0,1], doNorm=True,  n_bins=20, sepGenOpt='sepGenSig;sepGenBkg', doCuts=True, addData=False)
    #StackedHist(cfg.MC_samples,'NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=0.8', doCuts=True, addData=True)

    eff_tag  = 'tight'

    mu_eff     = AnaDict.read_pickle(f'files/{y}/eff_files/muon_eff_{y}_{eff_tag}.pkl')
    ele_eff    = AnaDict.read_pickle(f'files/{y}/eff_files/electron_eff_{y}_{eff_tag}.pkl')
    #
    mu_pt_bins = [float(bin_str.split(',')[0]) for bin_str in mu_eff['pt']] + [float(list(mu_eff['pt'].keys())[-1].split(',')[1])]
    ele_pt_bins = [float(bin_str.split(',')[0]) for bin_str in ele_eff['pt']] + [float(list(ele_eff['pt'].keys())[-1].split(',')[1])]
    mu_eta_bins = [float(bin_str.split(',')[0]) for bin_str in mu_eff['eta']] + [float(list(mu_eff['eta'].keys())[-1].split(',')[1])]
    ele_eta_bins = [float(bin_str.split(',')[0]) for bin_str in ele_eff['eta']] + [float(list(ele_eff['eta'].keys())[-1].split(',')[1])]

    
    
    
    StackedHist(cfg.MC_pow,'Lep_pt', xlabel='Elec_pt_pass_trigger_electron_pbt (GeV)',   bin_range=[0,700],      bins=ele_pt_bins,    add_cuts='MET_pt>=20;passSingleLepElec==1;pbt_elec==1',    doCuts=False, addData=True)
    #StackedHist(cfg.MC_samples,'Lep_eta', xlabel='Elec_eta_pass_trigger_electron_pbt', bin_range=[-2.6,2.6],   bins=ele_eta_bins,   add_cuts='MET_pt>=20;passSingleLepElec==1;pbt_elec==1',    doCuts=False, addData=True)
    #
    #
    StackedHist(cfg.MC_pow,'Lep_pt',   xlabel='Mu_pt_pass_trigger_muon_pbt (GeV)',  bin_range=[0,700],    bins=mu_pt_bins,   add_cuts='MET_pt>=20;passSingleLepMu==1;pbt_muon==1',  doCuts=False, addData=True)
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

