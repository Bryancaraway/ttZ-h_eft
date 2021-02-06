
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
from modules.AnaDict import AnaDict
from lib.fun_library import save_pdf
import config.ana_cff as cfg
###----=======-----###

#sepGenOpt ; 'sepGenSig','sepGenBkg','sepGenMatchedSig','sepGenMatchedBkg'

###----=======-----###

#jjec = 'ak8'
#jec_list = ['JESUp','JESDown','JERUp','JERDown']
processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','single_t','VJets','other']

@save_pdf('control_plots_NN.pdf')
#@save_pdf('doubleb_plots.pdf')
#@save_pdf('NN_compare.pdf')
#@save_pdf('high_eta_lep.pdf')
def main():
    for y in cfg.Years: 
        #for jec in jec_list:
        #Plotter.load_data(y, addBSF=False, tag=f'{jjec}{jec}') #tag='ak4JESUp'
        Plotter.load_data(y, samples=cfg.Sig_MC+cfg.Bkg_MC, addBSF=False, byprocess=True)
        ''' LOOK AT STACKED DATA VS MC '''
        #StackedHist(processes,    'Lep_pt', xlabel= r'lepton $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,   doCuts=False, add_cuts='MET_pt>20',   addData=True,  doShow=False)  
        ##StackedHist(processes,    'Lep_pt', xlabel= r'lepton $p_{T}$ (GeV) lep_eta>2', bin_range=[20,750],    n_bins=30, add_cuts='Lep_eta>2;MET_pt>20',   doCuts=False, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'lepton $\eta$',        bin_range=[-2.6,2.6],  n_bins=26,   doCuts=False, add_cuts='MET_pt>20', addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_pt', xlabel= r'Electron $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepElec==1;MET_pt>20', doCuts=False, addData=True,  doShow=False)  
        ##StackedHist(processes,    'Lep_pt', xlabel= r'Electron $p_{T}$ (GeV) lep_eta>2', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepElec==1;Lep_eta>2;MET_pt>20', doCuts=False, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'Electron $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepElec==1;MET_pt>20',  doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_pt', xlabel= r'Muon $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepMu==1;MET_pt>20', doCuts=False, addData=True,  doShow=False)  
        ##StackedHist(processes,    'Lep_pt', xlabel= r'Muon $p_{T}$ (GeV) lep_eta>2', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepMu==1;Lep_eta>2;MET_pt>20', doCuts=False, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'Muon $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepMu==1;MET_pt>20',  doCuts=False, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'MET_pt', xlabel=r'missing $e_{T}$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'nBottoms', xlabel='# of ak4 b-jets', bin_range=[-0.5,8.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],  doCuts=False, addData=True, doShow=False) 
        #StackedHist(processes,    'n_ak4jets', xlabel='# of ak4 jets',  bin_range=[-0.5,12.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5],  doCuts=False, addData=True, doShow=False) 
        #StackedHist(processes,    'n_ak8jets', xlabel='# of ak8 jets', bin_range=[-0.5,6.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5],  doCuts=False, addData=True, doShow=False) 
        ##
        #StackedHist(processes,    'jetpt_1', xlabel=r'leading jet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'jetpt_2', xlabel=r'next-to-leading jet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetpt_1', xlabel=r'leading bjet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetpt_2', xlabel=r'next-to-leading bjet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'fjetpt_1', xlabel=r'leading fatjet $p_{T} (AK8)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'jeteta_1', xlabel=r'leading jet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'jeteta_2', xlabel=r'next-to-leading jet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'bjeteta_1', xlabel=r'leading bjet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'bjeteta_2', xlabel=r'next-to-leading bjet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'fjeteta_1', xlabel=r'leading fatjet $\eta (AK8)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=False, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'jetbtag_1', xlabel=r'leading jet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'jetbtag_2', xlabel=r'next-to-leading jet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetbtag_1', xlabel=r'leading bjet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'bjetbtag_2', xlabel=r'next-to-leading bjet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        ##
        #StackedHist(processes,    'fjetsdm_1', xlabel=r'leading fatjet $M_{sd}$ (GeV)', bin_range=[50,200],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'fjetbbvl_1', xlabel=r'leading fatjet deepAK8MD bbvL', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'fjetwscore_1', xlabel=r'leading fatjet deepAK8MD WvsQCD', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,    'fjettscore_1', xlabel=r'leading fatjet deepAK8MD TvsQCD', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  
        #
        #StackedHist(processes,    'Zh_pt', bin_range=[200,500],  n_bins=20, add_cuts='Zh_bbvLscore>0.8',  doCuts=True, addData=True, doShow=False)  
        ###   #
        #StackedHist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='Zh_bbvLscore>0.8',  doCuts=True, addData=True, doShow=False)  
        #
        #StackedHist(processes,    'Lep_pt', xlabel= r'lepton $p_{T}$ (GeV) lep_eta>2, tight baseline', bin_range=[20,750],    n_bins=30, add_cuts='Lep_eta>2;MET_pt>20',   doCuts=True, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'lepton $\eta$, tight baseline',        bin_range=[-2.6,2.6],  n_bins=26,   doCuts=True, add_cuts='MET_pt>20', addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_pt', xlabel= r'Electron $p_{T}$ (GeV) lep_eta>2. tight baseline', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepElec==1;Lep_eta>2;MET_pt>20', doCuts=True, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'Electron $\eta$ tight baseline',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepElec==1;MET_pt>20',  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,    'Lep_pt', xlabel= r'Muon $p_{T}$ (GeV) lep_eta>2, tight baseline', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepMu==1;Lep_eta>2;MET_pt>20', doCuts=True, addData=True,  doShow=False)  
        #StackedHist(processes,    'Lep_eta', xlabel=r'Muon $\eta$ tight baseline',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepMu==1;MET_pt>20',  doCuts=True, addData=True, doShow=False)  
        # ----
        #StackedHist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLscore (no extra cuts)', bin_range=[0,1],  n_bins=20,  doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLscore Fake?', bin_range=[0,1],  n_bins=20, add_cuts='Zh_M>65;Zh_M<90;Zh_closeb_invM>150;Zh_closeb_invM<190',  doCuts=False, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bbvLscore', xlabel='Zh_bbvLscore Fake? (>2 b-jets outZH)', bin_range=[0,1],  n_bins=20, add_cuts='Zh_M>65;Zh_M<90;Zh_closeb_invM>150;Zh_closeb_invM<190;n_b_outZh>=2',  doCuts=False, addData=True, doShow=False)  
        
        StackedHist(processes,'withdak8md_NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='withdak8md_NN<=0.7',  doCuts=True, addData=True, doShow=False)  
        StackedHist(processes,'noak8md_NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='noak8md_NN<=0.7',     doCuts=True, addData=True, doShow=False)  
        StackedHist(processes,'withbbvl_NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', add_d_cuts='withbbvl_NN<=0.7',     doCuts=True, addData=True, doShow=False)  
        ###StackedHist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='NN>0.80;Zh_pt>300', doCuts=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withdak8md_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withdak8md_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (noak8md_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='noak8md_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
        #Hist(processes,'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) (withbbvl_NN>0.80,ZhpT>300)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='withbbvl_NN>0.80;Zh_pt>300', doCuts=True,  doNorm=False, doLog=True, addData=False, doShow=False)  
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
        #StackedHist(processes,'Zh_worstb_sj', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bestb_sj', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_bbscore_sj', bin_range=[0,2],  n_bins=20,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_Tscore', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
        #StackedHist(processes,'Zh_Wscore', bin_range=[0,1],  n_bins=10,  doCuts=True, addData=True, doShow=False)  
    
    #
    return 1



if __name__ == '__main__':
    main()
