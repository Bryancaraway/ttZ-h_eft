from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import numpy as np
from modules.AnaDict import AnaDict
from lib.fun_library import save_pdf, getFakebbvlCuts, getFakebbvlWeights, getWeightsWithEFT
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
#@save_pdf('control_plots_tight.pdf')
#@save_pdf('EFTLO_vs_NLO.pdf')
#@save_pdf('EFTLO_vs_NLO_templates_ttbb.pdf')
#@save_pdf('EFTLO_vs_NLO_mass_ttbb.pdf')
#@save_pdf('single_top_mcstats_v2.pdf')
#@save_pdf('control_plots_genzhpt.pdf')
#@save_pdf('sigbkg_purity.pdf')
#@save_pdf('Zh_info.pdf')
#@save_pdf('ttbb_lheht.pdf')
#@save_pdf('hl2_outputs.pdf')
#@save_pdf('ttzbb_efficiency.pdf')
#@save_pdf('NN_compare.pdf')
#@save_pdf('mass_compare.pdf')
#@save_pdf('fakebbvlsf_CR.pdf')
#@save_pdf('fakebbvlsf_CR_Andrew.pdf')
#@save_pdf("met_withqcd.pdf")
#@save_pdf("qcd_study_by_year_alt.pdf")
#@save_pdf("nn_comparison.pdf")
#@save_pdf("ttx_contamination.pdf")
def main():
    #for y in cfg.Years: 
    for y in ['2018']: 
        print(y)
        #for jec in jec_list:
        #Plotter.load_data(y, addBSF=False, tag=f'{jjec}{jec}') #tag='ak4JESUp'
        Plotter.load_data(y, samples=cfg.Sig_MC+cfg.Bkg_MC+['QCD'], addBSF=False, byprocess=True)
        #Plotter.load_data(y, samples=['single_t'], addBSF=False, byprocess=True)
        #Plotter.load_data(y, samples=cfg.Sig_MC+cfg.Bkg_MC+["QCD"]+['TTZ_EFT','TTH_EFT','TTbb_EFT'], addBSF=False, byprocess=True)
        #Plotter.load_data(y, samples=cfg.Sig_MC+['ttbb']+['TTZ_EFT','TTH_EFT','TTbb_EFT'], addBSF=False, byprocess=True)
        ''' LOOK AT STACKED DATA VS MC '''
        # AN, Paper Draft
        make_control_plots_tight()
        #make_control_plots_anastrat()
        #make_met_withqcd()
        #make_NN_compare()
        #make_qcd_study_by_year_alt()
        # = ==== = # 
        # extra
        #make_plots_ttX_cont()
        #make_plots_EFTLO_vs_NLO(y)
        #make_plots_singlet_stats()
        #make_genZhpt_plots()
        #make_sigbkg_purity()
        
    return 1
        

def make_control_plots_tight():
    # --- control plots tight
    Hist(processes,    cfg.nn, xlabel= r'NN with bbvl match', bin_range=[0,1],    n_bins=10,      doNorm=False, doCuts=True, add_cuts='bbvl_genmatch==1', doShow=True)  
    Hist(processes,    cfg.nn, xlabel= r'NN with no bbvl match', bin_range=[0,1],    n_bins=10,   doNorm=False, doCuts=True, add_cuts='bbvl_genmatch==0', doShow=True)  
    #StackedHist(processes,    'HT', xlabel= r'HT', bin_range=[0,2000],    n_bins=35,   doCuts=True,  addData=True, doShow=True)  
    exit()
    StackedHist(processes,    'PV_npvsGood', xlabel= r'nPVs', bin_range=[0,70],    n_bins=35,   doCuts=True,  addData=True, doShow=False)  
    #
    StackedHist(processes,    'Lep_pt', xlabel= r'lepton $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,   doCuts=True,  addData=True, doShow=False)  
    StackedHist(processes,    'Lep_eta', xlabel=r'lepton $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'Lep_pt', xlabel= r'Electron $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepElec==1', doCuts=True, addData=True,  doShow=False)  
    StackedHist(processes,    'Lep_eta', xlabel=r'Electron $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepElec==1',  doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'Lep_pt', xlabel= r'Muon $p_{T}$ (GeV)', bin_range=[20,750],    n_bins=30,  add_cuts='passSingleLepMu==1', doCuts=True, addData=True,  doShow=False)  
    StackedHist(processes,    'Lep_eta', xlabel=r'Muon $\eta$',        bin_range=[-2.6,2.6],  n_bins=26, add_cuts='passSingleLepMu==1',  doCuts=True, addData=True, doShow=False)  
    #
    StackedHist(processes,    'MET_pt', xlabel=r'missing $e_{T}$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'MET_phi', xlabel=r'missing $\phi$ (GeV)', bin_range=[-3.15,3.15],  n_bins=20,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'nBottoms', xlabel='# of ak4 b-jets', bin_range=[-0.5,8.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],  doCuts=True, addData=True, doShow=False) 
    StackedHist(processes,    'n_ak4jets', xlabel='# of ak4 jets',  bin_range=[-0.5,12.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5],  doCuts=False, addData=True, doShow=False) 
    StackedHist(processes,    'n_ak8jets', xlabel='# of ak8 jets', bin_range=[-0.5,6.5],  bins=[-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5],  doCuts=True, addData=True, doShow=False) 
    #
    StackedHist(processes,    'jetpt_1', xlabel=r'leading jet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'jetpt_2', xlabel=r'next-to-leading jet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'bjetpt_1', xlabel=r'leading bjet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'bjetpt_2', xlabel=r'next-to-leading bjet $p_{T} (AK4)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'fjetpt_1', xlabel=r'leading fatjet $p_{T} (AK8)$ (GeV)', bin_range=[0,500],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    #
    StackedHist(processes,    'jeteta_1', xlabel=r'leading jet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'jeteta_2', xlabel=r'next-to-leading jet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'bjeteta_1', xlabel=r'leading bjet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'bjeteta_2', xlabel=r'next-to-leading bjet $\eta (AK4)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'fjeteta_1', xlabel=r'leading fatjet $\eta (AK8)$', bin_range=[-2.6,2.6],  n_bins=26,     doCuts=True, addData=True, doShow=False)  
    #
    StackedHist(processes,    'jetbtag_1', xlabel=r'leading jet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'jetbtag_2', xlabel=r'next-to-leading jet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'bjetbtag_1', xlabel=r'leading bjet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    StackedHist(processes,    'bjetbtag_2', xlabel=r'next-to-leading bjet deepCSV (AK4)', bin_range=[0,1],  n_bins=25,     doCuts=True, addData=True, doShow=False)  
    #
    StackedHist(processes,    'fjetsdm_1', xlabel=r'leading fatjet $m_{sd}$ (GeV)', bin_range=[50,200],  n_bins=25,     doCuts=True, addData=True, doShow=False)  

def make_control_plots_anastrat():
    # ---- ana strat
    StackedHist(processes,    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV]', bins=[200,300,450,600],  doCuts=True, doLog=True, addData=True, doShow=False)  
    StackedHist(processes,'newgenm_NN', xlabel='DNN score', bin_range=[0,1],  n_bins=20,  add_d_cuts='newgenm_NN<=1.8',     doCuts=True, addData=True, doShow=False)  
    Hist(['ttZ','ttH','TTBar','tt_B'],'newgenm_NN', xlabel='DNN score', bin_range=[0,1],  n_bins=10, doNorm=True, doLog=True, doCuts=True, addData=False, doShow=False)  
    StackedHist(processes,'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV]', bin_range=[50,200],  bins=cfg.sdm_bins, doCuts=True,  doLog=True, addData=True, doShow=False)  
       #
    Hist(['ttZ','ttH','TTBar','tt_B'],'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV] (Z/H $p_{\mathrm{T}}>300$ (GeV);DNN$>0.8$)', bin_range=[50,200],  bins=[50,75,90,105,120,140,200], add_cuts='Zh_pt>300;newgenm_NN>0.80', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  
    Hist(['ttZ','ttH','TTBar','tt_B'],'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV] (Z/H $p_{\mathrm{T}}>450$ (GeV);DNN$>0.8$)', bin_range=[50,200],  bins=[50,75,90,105,120,140,200], add_cuts='Zh_pt>450;newgenm_NN>0.80', doCuts=True,  doNorm=True, doBinRatio=False, doLog=False, addData=False, doShow=False)  
     
def make_NN_compare():
    #StackedHist(processes,'newgenm_NN', xlabel='NN: new 50 vars (Zh_pt>300, Z/H M)', bins=[0. ,0.091, 0.343 ,0.459, 0.612, 0.753, 1.], add_cuts="Zh_pt>300;Zh_M>75;Zh_M<145", doLog=True, addData=True, doCuts=True, doShow=False)  
    #StackedHist(processes,'newreduced1p0_NN', xlabel='NN: new 28 vars (Zh_pt>300, Z/H M)', bins=[0. ,0.107, 0.306 ,0.385, 0.518, 0.662, 1.], add_cuts="Zh_pt>300;Zh_M>75;Zh_M<145", doLog=True, addData=True, doCuts=True, doShow=False)  
    Hist(processes,'newgenm_NN',       xlabel='NN: new 50 vars (Zh_pt>300, Z/H M)', bins=[0. ,0.121, 0.404, 0.540, 0.703, 0.832, 1.], add_cuts="Zh_pt>300;Zh_M>75;Zh_M<145", doLog=True, doNorm=False, addSoB=True, doCuts=True, doShow=False)  
    Hist(processes,'newreduced1p0_NN', xlabel='NN: new 28 vars (Zh_pt>300, Z/H M)', bins=[0. , 0.166, 0.430, 0.522, 0.669, 0.796, 1.], add_cuts="Zh_pt>300;Zh_M>75;Zh_M<145", doLog=True, doNorm=False, addSoB=True, doCuts=True, doShow=False)  
    #
#@save_pdf("met_withqcd.pdf")
def make_met_withqcd():
    StackedHist(processes+['QCD'],'MET_pt', xlabel=r'missing $e_{T}$ (GeV), with j/psi lepton cut and sip3d < 4', bin_range=[0,500],  n_bins=25,  doCuts=True, addData=True, doShow=False)  
#@save_pdf("qcd_study_by_year_alt.pdf") 
def make_qcd_study_by_year_alt():
    ### QCD study
    StackedHist(processes+['QCD'],cfg.nn, xlabel=r'MVA score with j/psi lepton cut and sip3d < 4', bin_range=[0,1],  n_bins=20, add_d_cuts=f'{cfg.nn}<=0.8',  doCuts=True, addData=True, doShow=False)  
    #Hist(['QCD','ttZ','ttH'],'min_sel_soft_mu_invm',xlabel='min_sel_soft_mu_invm', bin_range=[0,100],  n_bins=100, sepGenOpt='sepGenMatchedSig;+',   add_cuts='passSingleLepMu==1;NN<=1.80',    doCuts=True, doLog=True, doNorm=False, doShow=False)
    ##
    #Hist(['QCD','ttZ','ttH'],'Muon_sip3d',xlabel='Muon_sip3d', bin_range=[0,10],  n_bins=40, sepGenOpt='sepGenMatchedSig;+',  add_cuts='passSingleLepMu==1;NN<=1.80;passNotHadLep==1',    doCuts=True, doLog=True, doNorm=False, doShow=False)
    ##
    #Hist(['QCD','ttZ','ttH'],'min_sel_soft_elec_invm',xlabel='min_sel_soft_elec_invm', bin_range=[0,100],  n_bins=100, sepGenOpt='sepGenMatchedSig;+',   add_cuts='passSingleLepElec==1;NN<=1.80',    doCuts=True, doLog=True, doNorm=False, doShow=False)
    #Hist(['QCD','ttZ','ttH'],'Electron_sip3d',xlabel='Electron_sip3d', bin_range=[0,10],  n_bins=40, sepGenOpt='sepGenMatchedSig;+',  add_cuts='passSingleLepElec==1;NN<=1.80;Zh_pt>300;passNotHadLep==1',    doCuts=True, doLog=True, doNorm=False, doShow=False)
    #
    StackedHist(processes+['QCD'],'min_sel_soft_elec_invm', bin_range=[0,100], n_bins=100, add_cuts='passSingleLepElec==1', doCuts=False, doLog=True, addData=True, doShow=False)
    StackedHist(processes+['QCD'],'min_sel_soft_mu_invm', xlabel='min_sel_soft_mu_invm', bin_range=[0,100], n_bins=100, add_cuts='passSingleLepMu==1', doCuts=False, doLog=True, addData=True, doShow=False)
    #Hist(['QCD','ttZ','ttH'],'MET_pt', xlabel=r'missing $e_{T}$ (GeV), with j/psi lepton cut', bin_range=[0,500],  n_bins=25,  add_cuts='passNotHadLep==1', sepGenOpt='sepGenMatchedSig;+',   doCuts=True, doShow=False)  
    #Hist(['QCD','ttZ','ttH'],'MET_pt', xlabel=r'missing $e_{T}$ (GeV), with j/psi lepton cut, Zhpt>300', bin_range=[0,500],  n_bins=25,  add_cuts='passNotHadLep==1;Zh_pt>300', sepGenOpt='sepGenMatchedSig;+',   doCuts=True, doShow=False)  
    
#

def make_plots_ttX_cont():
    Hist(['ttX'], cfg.nn, xlabel=r'NN score', bin_range =[0,1],  n_bins=20, sepGenOpt='sepBySample',   doCuts=True, doShow=False)  


def make_plots_EFTLO_vs_NLO(y):
    alt_weight = getWeightsWithEFT
    nn = cfg.nn
    #sig_processes = ['ttZ','ttH','tt_B']
    sig_processes = ['tt_B']
    common_kwargs =  {'alt_weight':alt_weight, 'sepGenOpt':'sepByEFT', 'doCuts':True, 'doLog':False, 'doNorm':True, 'addData':False, 'doShow':False}
    nn_bins = { 
        '2016': {
            '200':[0.0,0.21,0.59,0.70,0.78,0.90,1.00],
            '300':[0.0,0.32,0.67,0.78,0.86,0.92,1.00],
            '450':[0.0,0.40,0.72,0.80,0.87,0.92,1.00],
        },
        '2017': {
            '200':[0.0,0.18,0.53,0.65,0.80,0.88,1.00],
            '300':[0.0,0.23,0.62,0.72,0.83,0.91,1.00],
            '450':[0.0,0.34,0.71,0.79,0.87,0.92,1.00],
        },
        '2018': {
            '200':[0.0,0.17,0.51,0.63,0.78,0.88,1.00],
            '300':[0.0,0.21,0.58,0.70,0.82,0.90,1.00],
            '450':[0.0,0.31,0.69,0.78,0.87,0.92,1.00],
        },
    }
    pt_bins = [200,300,450, np.inf]
    mass_bins = [50,75,105,145,200]
    for sig in sig_processes:
        #Hist([sig],    'ttbb_genbb_invm',  xlabel='ttbb_genbb_invm ', bin_range=[0,200], n_bins=20, add_cuts=f'{nn}>0.0', **common_kwargs)
        #Hist([sig],    'ttbb_genbb_invm',  xlabel='ttbb_genbb_invm (NN > 0.8)', bin_range=[0,200], n_bins=20, add_cuts=f'{nn}>0.8', **common_kwargs)
        #Hist([sig],    'ttbb_genbb_invm',  xlabel='ttbb_genbb_invm (NN > 0.8,Zh_pt>300,Zh_M>75,Zh_M<105)', bin_range=[0,200], n_bins=20, add_cuts=f'{nn}>0.8;Zh_pt>300;Zh_M>75;Zh_M<105', **common_kwargs)
        #Hist([sig],    'ttbb_genbb_pt',    xlabel='ttbb_genbb_pt ', bin_range=[0,600], n_bins=20, add_cuts=f'{nn}>0.0', **common_kwargs)
        #Hist([sig],    'ttbb_genbb_pt',    xlabel='ttbb_genbb_pt (NN > 0.8)', bin_range=[0,600], n_bins=20, add_cuts=f'{nn}>0.8', **common_kwargs)
        #Hist([sig],    'ttbb_genbb_pt',    xlabel='ttbb_genbb_pt (NN > 0.8,Zh_pt>300,Zh_M>75,Zh_M<105)', bin_range=[0,600], n_bins=20, add_cuts=f'{nn}>0.8;Zh_pt>300;Zh_M>75;Zh_M<105', **common_kwargs)
        #Hist([sig],    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV]', bins=[200,300,450,600], **common_kwargs)
        #Hist([sig], nn, xlabel='DNN score', bin_range=[0,1],  n_bins=10, **common_kwargs)
        #Hist([sig],'Zh_M', xlabel=r'Z/H $m_{\mathrm{SD}}$ [GeV]', bin_range=[50,200],  bins=cfg.sdm_bins, **common_kwargs)  
        for i_bin in range(1,len(pt_bins)):
            for j_bin in range(1,len(mass_bins)):
                #Hist([sig], nn, xlabel=f'DNN score (Z/H pT {i_bin}, Z/H mass {j_bin})', bins=nn_bins[y][str(pt_bins[i_bin-1])], 
                Hist([sig], nn, xlabel=f'DNN score (Z/H pT {i_bin}, Z/H mass {j_bin})', bins=[0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0], 
                     add_cuts=f'Zh_pt>{pt_bins[i_bin-1]};Zh_pt<{pt_bins[i_bin]};Zh_M>{mass_bins[j_bin-1]};Zh_M<{mass_bins[j_bin]}', **common_kwargs)
    

def make_plots_singlet_stats():
    Hist(['single_t'], cfg.nn, xlabel=r'NN score (Z/H cand pT > 450)', bin_range =[0,1],  n_bins=10, add_cuts='Zh_pt>450', sepGenOpt='sepBySample',  doNorm=False,  doCuts=True, doShow=False)  
    Hist(['single_t'], cfg.nn, xlabel=r'NN score (Z/H cand pT > 450, Z mass)', bins =[0.0,0.34,0.71,0.79,0.87,0.92,1.0],  add_cuts='Zh_pt>450;Zh_M>75;Zh_M<105', sepGenOpt='sepBySample',  doNorm=False,  doCuts=True, doShow=False)  
    Hist(['single_t'], cfg.nn, xlabel=r'NN score (n_w_leps > 0, Z/H pt > 300)', bin_range =[0,1],  n_bins=10, add_cuts='Zh_pt>300;n_w_leps>0',   doNorm=False,  doCuts=True, doShow=False)  
    Hist(['single_t'], cfg.nn, xlabel=r'NN score (n_w_leps = 0, Z/H pt > 300)', bin_range =[0,1],  n_bins=10, add_cuts='Zh_pt>300;n_w_leps==0',   doNorm=False,  doCuts=True, doShow=False)  

def make_genZhpt_plots():
    Hist(['ttH','ttZ'],    'genZHpt', xlabel=r'GEN Z/H $p_{T}$ [GeV]', n_bins=30, bin_range=[0,600], doNorm=True,  doCuts=False, doLog=True, doShow=False)      

def make_sigbkg_purity():
    #Hist(['ttH','ttZ'],    'genZHpt', xlabel=r'GEN Z/H $p_{T}$ [GeV]', bins=[0,200,300,450,600], doNorm=True,  doCuts=True, doLog=False, doShow=False)      
    Hist(['TTBar'],        'tt_C', xlabel=r'is tt_C', bins=[-0.5,0.5,1.5], doNorm=True, doCuts=True,  doLog=False, doShow=False)      
    Hist(['tt_B'],        'tt_2b', xlabel=r'is tt_2b', bins=[-0.5,0.5,1.5], doNorm=True, doCuts=True, doLog=False, doShow=False)      
    #
    #Hist(['ttH','ttZ'],    'genZHpt', xlabel=r'GEN Z/H $p_{T}$ [GeV] (NN >0.8)', bins=[0,200,300,450,600], add_cuts=f'{cfg.nn}> 0.8', doNorm=True,  doCuts=True, doLog=False, doShow=False)      
    Hist(['TTBar'],        'tt_C',    xlabel=r'is tt_C (NN >0.8)',               bins=[-0.5,0.5,1.5],      add_cuts=f'{cfg.nn}> 0.8', doNorm=True, doCuts=True,  doLog=False, doShow=False)      
    Hist(['tt_B'],        'tt_2b',    xlabel=r'is tt_2b (NN >0.8)',              bins=[-0.5,0.5,1.5],      add_cuts=f'{cfg.nn}> 0.8', doNorm=True, doCuts=True, doLog=False, doShow=False)      
    

###########
#StackedHist(processes,    'fjetbbvl_1', xlabel=r'leading fatjet deepAK8MD bbvL', bin_range=[0,1],  n_bins=25,     doCuts=False, addData=True, doShow=False)  

###   #
#StackedHist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV)', bin_range=[50,200],  bins=[50,80,105,145,200], add_cuts='Zh_bbvLscore>0.8',  doCuts=True, addData=True, doShow=False)  
#StackedHist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV)', bin_range=[50,200], n_bins=50, add_cuts='n_b_outZh>=2;Zh_bbvLscore<0.4', doLog=False,  doCuts=False, addData=True, doShow=True)  

#Hist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ [GeV] (Zhpt>300,NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='NN>0.8;Zh_pt>300',  doCuts=True, addSoB=True, doShow=False)  
#Hist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ [GeV] (Zhpt>450,NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='NN>0.8;Zh_pt>450',  doCuts=True, addSoB=True, doShow=False)  
#Hist(['ttZ','tt_B','TTBar', 'single_t', 'ttX', 'VJets'],    'Zh_M', xlabel=r'Z $m_{sd}$ [GeV] (Zhpt>300,NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='NN>0.8;Zh_pt>300',  doCuts=True, addSoB=True, doShow=False)  
#Hist(['ttH','tt_B','TTBar', 'single_t', 'ttX', 'VJets'],    'Zh_M', xlabel=r'H $m_{sd}$ [GeV] (Zhpt>300,NN>0.8)', bin_range=[50,200],  bins=30, add_cuts='NN>0.8;Zh_pt>300',  doCuts=True, addSoB=True, doShow=False)  
##
###------- mass validation
#StackedHist(processes,    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) baseline, pt>300', bin_range=[50,200],      n_bins=30, add_cuts='Zh_pt>300',  doCuts=True, addData=True, doShow=False)  
#StackedHist(processes,    'Zh_M_nom', xlabel=r"Z/H $m_{sd}$ 'bad' (GeV) baseline, pt>300", bin_range=[50,200],  n_bins=30, add_cuts='Zh_pt>300',  doCuts=True, addData=True, doShow=False)  
#StackedHist(processes,    'Zh_M_alt', xlabel=r"Z/H $m_{sd}$ sjcorr (GeV) baseline, pt>300", bin_range=[50,200],  n_bins=30, add_cuts='Zh_pt>300',  doCuts=True, addData=True, doShow=False)  
#StackedHist(processes,    'fjetsdm_1', xlabel=r'fj1 $m_{sd}$ (GeV) wmd>0.5, pt>300', bin_range=[50,200],      n_bins=30, add_cuts='fjetwmdscore_1>0.5;fjetpt_1>300',  doCuts=False, addData=True, doShow=False)  
#StackedHist(processes,    'fjetsdmnom_1', xlabel=r"fj1 $m_{sd}$ 'bad' (GeV) wmd>0.5, pt>300", bin_range=[50,200],  n_bins=30, add_cuts='fjetwmdscore_1>0.5;fjetpt_1>300',  doCuts=False, addData=True, doShow=False)  
#StackedHist(processes,    'fjetsdmalt_1', xlabel=r"fj1 $m_{sd}$ sjcorr (GeV) wmd>0.5, pt>300", bin_range=[50,200],  n_bins=30, add_cuts='fjetwmdscore_1>0.5;fjetpt_1>300',  doCuts=False, addData=True, doShow=False)  
#
#Hist(['ttZ','ttH','tt_B','TTBar'],    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) baseline, pt 300-450', bin_range=[50,200],      n_bins=30, add_cuts='Zh_pt>300;Zh_pt<450',  doCuts=True, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'Zh_M_nom', xlabel=r"Z/H $m_{sd}$ 'bad' (GeV) baseline, pt 300-450", bin_range=[50,200],  n_bins=30, add_cuts='Zh_pt>300;Zh_pt<450',  doCuts=True, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'Zh_M_alt', xlabel=r"Z/H $m_{sd}$ sjcorr (GeV) baseline, pt 300-450", bin_range=[50,200],  n_bins=30, add_cuts='Zh_pt>300;Zh_pt<450',  doCuts=True, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'Zh_M', xlabel=r'Z/H $m_{sd}$ (GeV) baseline, pt 300-450, NN>0.8', bin_range=[50,200],      n_bins=30, add_cuts='Zh_pt>300;Zh_pt<450;NN>0.8',  doCuts=True, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'Zh_M_nom', xlabel=r"Z/H $m_{sd}$ 'bad' (GeV) baseline, pt 300-450, NN>0.8", bin_range=[50,200],  n_bins=30, add_cuts='Zh_pt>300;NN>0.8;Zh_pt<450',  doCuts=True, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'Zh_M_alt', xlabel=r"Z/H $m_{sd}$ sjcorr (GeV) baseline, pt 300-450, NN>0.8", bin_range=[50,200],  n_bins=30, add_cuts='Zh_pt>300;NN>0.8;Zh_pt<450',  doCuts=True, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'fjetsdm_1', xlabel=r'fj1 $m_{sd}$ (GeV) wmd>0.5, pt>300', bin_range=[50,200],      n_bins=30, add_cuts='fjetwmdscore_1>0.5;fjetpt_1>300',  doCuts=False, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'fjetsdmnom_1', xlabel=r"fj1 $m_{sd}$ 'bad' (GeV) wmd>0.5, pt>300", bin_range=[50,200],  n_bins=30, add_cuts='fjetwmdscore_1>0.5;fjetpt_1>300',  doCuts=False, addData=True, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],    'fjetsdmalt_1', xlabel=r"fj1 $m_{sd}$ sjcorr (GeV) wmd>0.5, pt>300", bin_range=[50,200],  n_bins=30, add_cuts='fjetwmdscore_1>0.5;fjetpt_1>300',  doCuts=False, addData=True, doShow=False)  

###------- end of mass validation
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
#Hist(['ttZ','ttH'],'n_tt_leps_notau', bins=[-0.5,0.5,1.5,2.5], add_cuts='Zllnunu==0;Hnonbb==0', doNorm=True, doCuts=True, doShow=True)  

#----------- N lep fraction -------- #

#Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps', xlabel='n_tt_leps', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, doNorm=True, doShow=False)
#Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps', xlabel='n_tt_leps (NN>0.8)', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, add_cuts='NN>0.8', doNorm=True, doShow=False)
#
#Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps_notau', xlabel='n_tt_leps no taus', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, doNorm=True, doShow=False)
#Hist(['ttZ','ttH', 'TTBar'], 'n_tt_leps_notau', xlabel='n_tt_leps no taus (NN>0.8)', bin_range=[-0.5,3.5], n_bins=4, doCuts=True, add_cuts='NN>0.8', doNorm=True, doShow=False)

# --------------------------- #


#Hist(['ttZ','ttH'],'NN',xlabel='NN', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80', sepGenOpt='sepGenSig', doCuts=True, doLog=True, doNorm=False, doShow=True)  
#Hist(['ttZ','ttH'],'withbbvl_NN',xlabel='NN (Zh_pt>300)', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80;Zh_pt>300', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=False)  
#Hist(['ttZ','ttH','tt_B','TTBar'],'NN',xlabel=r'DNN score (${p}_{\mathrm{T}}^{\mathrm{Z/H}}>300$ GeV)', bin_range=[0,1],  n_bins=20, add_cuts='NN<=1.80;Zh_pt>300', sepGenOpt='sepGenMatchedSig;+', doCuts=True, doLog=True, doNorm=True, doShow=False)  
#Hist(['ttZ','ttH'],'Zh_bbvLscore',xlabel='bbvL', bin_range=[.8,1],  n_bins=5, add_cuts='NN<=1.80', sepGenOpt='sepGenMatchedSig;++', doCuts=True, doNorm=False, doShow=True)  
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
#StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$', bin_range=[-3.15,3.15], n_bins=20,  doCuts=True, doLog=False, addData=True, doShow=False)  
#StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$ (Z/H $\eta < 0$)', bin_range=[-3.15,3.15], n_bins=20, add_cuts='Zh_eta<0', doCuts=True, doLog=False, addData=True, doShow=False)  
#StackedHist(processes,    'Zh_phi', xlabel=r'Z/H $\phi$ (Z/H $\eta > 0$)', bin_range=[-3.15,3.15], n_bins=20, add_cuts='Zh_eta>0', doCuts=True, doLog=False, addData=True, doShow=False)  
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

# HT seperation for Ken

#Hist(['tt_B'],    'LHE_HTIncoming', xlabel='LHE_HTIncoming (inc)', n_bins=50, bin_range=[0, 1500], doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
#Hist(['tt_B'],    'LHE_HTIncoming', xlabel='LHE_HTIncoming (300 < ZHpt < 450)', n_bins=50, bin_range=[0, 1500], add_cuts='Zh_pt>300;Zh_pt<450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
#Hist(['tt_B'],    'LHE_HTIncoming', xlabel='LHE_HTIncoming (ZHpt > 450)', n_bins=50, bin_range=[0, 1500],       add_cuts='Zh_pt>450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
#
#Hist(['tt_B'],    'LHE_HT', xlabel='LHE_HT (inc)', n_bins=50, bin_range=[0, 1500], doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
#Hist(['tt_B'],    'LHE_HT', xlabel='LHE_HT (300 < ZHpt < 450)', n_bins=50, bin_range=[0, 1500], add_cuts='Zh_pt>300;Zh_pt<450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  
#Hist(['tt_B'],    'LHE_HT', xlabel='LHE_HT (ZHpt > 450)', n_bins=50, bin_range=[0, 1500],       add_cuts='Zh_pt>450',doCuts=True, doLog=True, doNorm=False, addData=False,  doShow=False)  

#Hist(['ttZ','ttH'],    'Zh_pt', xlabel=r'Z/H $p_{\mathrm{T}}$ [GeV]', bins=[200,300,450,600], sepGenOpt='sepGenMatchedSig;++', doCuts=True,  doLog=True, doNorm=False, doShow=False)  
#Hist(['ttZ','ttH'],    'matchedGen_ZHbb_tt', bins=[-0.5, 0.5, 1.5, 2.5], doCuts=True, sepGenOpt='sepGenSig', doLog=True, doNorm=False, doShow=False)  

# NN hidden layer # 2 output

#for i in range(64): # 64 perceptrons in hidden layer # 2
#    StackedHist(processes,f'NN_{i}', xlabel=f'hlayer#2, node {i}', bin_range=[-3,3],  n_bins=20,  doCuts=True, addData=True, doShow=False)              


if __name__ == '__main__':
    main()
