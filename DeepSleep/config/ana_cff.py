################
## Config for ##
# TTV analysis #
################
#
import subprocess as sb
# get current working directory according to git
def cdir():
    _cdir = sb.check_output('pwd', shell=True).decode().strip('\n')
    if 'condor' in _cdir: _wdir = './' # if working on condor, already in the right directory
    else:
        _wdir = sb.check_output('echo $(git rev-parse --show-toplevel)', shell=True).decode().strip('\n')+'/DeepSleep/'
    return _wdir, _cdir

_wdir, _cdir = cdir()
#
master_file_path  = _wdir+'files/'
dataDir           = _wdir+'data/'
pdfDir            = _wdir+'pdf/'
# Overhead #
import os
if   os.path.exists('/cms/data/store/user/ttxeft/') : # test to see if on kodiak
    file_path         = '/cms/data/store/user/ttxeft/Skim_nanoAOD/' # for kodiak
elif os.path.exists('/eos/uscms/') or 'condor' in _cdir: # test to see if on lpc will need to fix for condor on kodiak i think
    file_path        = 'root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/'
else: raise("Not on Kodiak or LPC, please manually input file_path in file: ./config/ana_cff.py")

tree_dir          = 'Training'
##
ZHptcut           = 200
Years             = ['2016','2017','2018']
#MC_samples        = ['TTZH', 'QCD', 'TTX', 'DY', 'WJets', 'TTBarHad', 'DiBoson', 'TriBoson', 'TTBarLep']#,'ZJets']

ttbar_samples     = ['TTBarHad_pow','TTBarSemi_pow','TTBarDi_pow','TTbbHad_pow','TTbbSemi_pow','TTbbDi_pow']
#['TTBarHad_pow', 'TTBarLep_pow','TT_bb_pow', 'TTBarHad', 'TTBarLep']
MC_pow            = ['TTZH', 'QCD', 'TTX', 'DY', 'WJets', 'DiBoson', 'TriBoson', 'TTBarHad_pow','TTBarSemi_pow','TTBarDi_pow','TTbbHad_pow','TTbbSemi_pow','TTbbDi_pow']
Sig_MC            = [#'TTZH', 
                     'TTZ_bb','TTZ','TTH']
## EFT ##
Sig_EFT_MC        = ['TTZ_EFT','TTH_EFT']
tt_eft_samples    = ['TTJets_EFT','TTBB_EFT']
#########
Bkg_MC            = [#'QCD', 
    'TTX', 'DY', 'WJets', 'DiBoson', 'TriBoson', 'TTBarHad_pow','TTBarSemi_pow','TTBarDi_pow','TTbbHad_pow','TTbbSemi_pow','TTbbDi_pow']
All_MC            = ['TTZH', 'TTZ_bb', 'TTZ', 'TTH', 'QCD', 'TTX', 'DY', 'WJets','DiBoson', 'TriBoson', 'TTBarHad_pow','TTBarSemi_pow','TTBarDi_pow','TTbbHad_pow','TTbbSemi_pow','TTbbDi_pow']
#['TTZH', 'TTZ_bb', 'QCD', 'TTX', 'DY', 'WJets', 'TTBarHad', 'TTBarHad_pow', 'DiBoson', 'TriBoson', 'TTBarLep','TTBarLep_pow', 'TT_bb_pow']
# Handle systematic sample docs
tt_sys_samples    = ['TTBarHad_pow_erdOn','TTBarHad_pow_UEUp','TTBarHad_pow_UEDown','TTBarHad_pow_hdampUp','TTBarHad_pow_hdampDown',
                     'TTBarSemi_pow_erdOn','TTBarSemi_pow_UEUp','TTBarSemi_pow_UEDown','TTBarSemi_pow_hdampUp','TTBarSemi_pow_hdampDown',
                     'TTBarDi_pow_erdOn','TTBarDi_pow_UEUp','TTBarDi_pow_UEDown','TTBarDi_pow_hdampUp','TTBarDi_pow_hdampDown',
                     'TTbbHad_pow_hdampUp','TTbbHad_pow_hdampDown',
                     'TTbbSemi_pow_hdampUp','TTbbSemi_pow_hdampDown',
                     'TTbbDi_pow_hdampUp','TTbbDi_pow_hdampDown']
tt_bb             = ['TTbbHad_pow','TTbbSemi_pow','TTbbDi_pow']
tt_bb_sys         = [  'TTbbHad_pow_hdampUp','TTbbHad_pow_hdampDown',
                       'TTbbSemi_pow_hdampUp','TTbbSemi_pow_hdampDown',
                       'TTbbDi_pow_hdampUp','TTbbDi_pow_hdampDown']
#
jec_variations    = [jtype+jec for jec in ['JESUp','JESDown','JERUp','JERDown'] for jtype in ['ak4','ak8']]
sig_sys_samples   = [sig+'_'+jec for sig in Sig_MC for jec in jec_variations]
bkg_sys_samples   = [bkg+'_'+jec for bkg in Bkg_MC for jec in jec_variations] + tt_sys_samples
all_sys_samples   = sig_sys_samples + bkg_sys_samples
#
Data_samples      = ['EleData','MuData']
Lumi              = {'2016': 35.917149,
                     '2017': 41.525338,
                     '2018': 59.72444,
                     '2018preHEM' : 21.1,
                     '2018postHEM': 38.6,
                     'run2': 137.166648,
                     'Total': 137.166648
                  } 
##
##############
##### TTZ, Z to bb CONFIG #####
ZHbbFitMinJets = 4
ZHbbFitMaxJets = 100
ZHbb_btagWP    = {'2016': 0.6321, # Med for 2016
                  '2017': 0.4941, # Med for 2017
                  '2018': 0.4148  # Med for 2018
                  }
# ttZ/H->bb SM x-section
ZHbbXsec = {'ttZbb': .1157,
            'ttHbb': .2934 }
ZHbbtotXsec = ZHbbXsec['ttZbb'] + ZHbbXsec['ttHbb']
# ttZ/H->bb MC count 2017
n_ZHbbMC_dict      = {'ttZbb': 163876,
                      'ttHbb': 5698653 }
n_ZHbbMC           = n_ZHbbMC_dict['ttZbb'] + n_ZHbbMC_dict['ttHbb']
#
hlt_path = {
    'muon'    :{ '2016': (lambda x : ((x['HLT_IsoMu24']) | 
                                      (x['HLT_IsoTkMu24']) | 
                                      (x['HLT_Mu50']) | 
                                      (x['HLT_TkMu50']))),

                 '2017': (lambda x : ((x['HLT_IsoMu27']) | 
                                      (x['HLT_Mu50']) | 
                                      (x['HLT_OldMu100']) | 
                                      (x['HLT_TkMu100']))),

                 '2018': (lambda x : ((x['HLT_IsoMu24']) | 
                                      (x['HLT_Mu50']) | 
                                      (x['HLT_OldMu100']) | 
                                      (x['HLT_TkMu100']))),
             },
    'electron':{ '2016': (lambda x : ((x['HLT_Ele27_WPTight_Gsf']) | 
                                      (x['HLT_Photon175']) | 
                                      (x['HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165']) | 
                                      (x['HLT_Ele115_CaloIdVT_GsfTrkIdT']) | 
                                      (x['HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50']))),

                 '2017': (lambda x : ((x['HLT_Ele32_WPTight_Gsf_L1DoubleEG']) | 
                                      (x['HLT_Ele35_WPTight_Gsf']) | 
                                      (x['HLT_Photon200']) | 
                                      (x['HLT_Ele115_CaloIdVT_GsfTrkIdT']) | 
                                      (x['HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165']) #|
                                      #(x['HLT_Ele28_eta2p1_WPTight_Gsf_HT150']) # trying this out
                                  )),

                 '2018': (lambda x : ((x['HLT_Ele32_WPTight_Gsf']) | 
                                      (x['HLT_Ele115_CaloIdVT_GsfTrkIdT']) | 
                                      (x['HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165']) | 
                                      (x['HLT_Photon200']) #|
                                      #(x['HLT_Ele28_eta2p1_WPTight_Gsf_HT150']) # trying this out
                                  ))
             }
}
###################
# Input Variables #
#LC = '_drLeptonCleaned'
LC = ''
#
lep_sel_vars = {'muon'    : ['Muon_pt','Muon_eta','Muon_phi','Muon_mass',
                             'Muon_miniPFRelIso_all','Muon_mediumId'],#'Muon_FlagId'],
                'electron': ['Electron_pt','Electron_eta','Electron_phi','Electron_mass',
                             'Electron_miniPFRelIso_all', 'Electron_cutBasedNoIso']}

lep_sel =      {'muon': (lambda x: ((x['Muon_pt'] > 30)        & 
                                    (abs(x['Muon_eta']) < 2.4) &
                                    #(x['Muon_FlagId'] >= 1)    & 
                                    (x['Muon_mediumId'] >= 1)    & 
                                    (x['Muon_miniPFRelIso_all'] < 0.2) )),

                'electron': {'2016': (lambda x : ((x['Electron_pt'] > 30) & (abs(x['Electron_eta']) < 2.5) & 
                                                  (x['Electron_cutBasedNoIso'] >= 4) & (x['Electron_miniPFRelIso_all'] < 0.1))),
                             '2017': (lambda x : ((x['Electron_pt'] > 35) & (abs(x['Electron_eta']) < 2.5) & 
                                                  (x['Electron_cutBasedNoIso'] >= 4) & (x['Electron_miniPFRelIso_all'] < 0.1))),
                             '2018': (lambda x : ((x['Electron_pt'] > 35) & (abs(x['Electron_eta']) < 2.5) & 
                                                  (x['Electron_cutBasedNoIso'] >= 4) & (x['Electron_miniPFRelIso_all'] < 0.1))),
                             }
                }
# JMS (jet mass scale) and JMR (jet mass resolution) for softdrop Mass
ak8_sys_vars = ['FatJet_pt'+LC, 'FatJet_eta'+LC, 'FatJet_phi'+LC, 'FatJet_mass'+LC, 
                'FatJet_msoftdrop'+LC, 'FatJet_rawFactor'+LC, 'FatJet_area'+LC,
                'GenJetAK8_pt', 'GenJetAK8_eta', 'GenJetAK8_phi', 'GenJetAK8_mass',
                #'SubGenJetAK8_pt', 'SubGenJetAK8_eta', 'SubGenJetAK8_phi', 'SubGenJetAK8_mass',
                #'SubJet_pt', 'SubJet_eta', 'SubJet_phi', 'SubJet_mass', 'SubJet_rawFactor',
                #'FatJet_subJetIdx1'+LC,'FatJet_subJetIdx2'+LC
            ]

ak8_softdropM_info = {'jms':{'value':0.999,
                             'up'   :0.999+0.004,
                             'down' :0.999-0.004},
                      'jmr':{'value':1.079,
                             'up'   :1.079+0.105,
                             'down' :1.079-0.105}
                  }
#ana_sf_dir    = '/cms/data/store/user/ttxeft/Ana_sf_files'

#why...
BC_btag_sf = {'2016': {'values': [1.0, 1.0, 0.94241418, 0.98314421, 1.05133896, 1.33768073, 0.54615971, 0.],                                     
                       'err':    [0, 0, 0.00209416, 0.00650319, 0.02603249, 0.13512616, 0.38619323, 0]},
              '2017': {'values': [1.0, 1.0, 0.97374629, 1.05901002, 1.22053029, 1.31213915, 1.38618462, 0.],
                       'err':    [0, 0, 0.00216191, 0.0071413, 0.03027763, 0.1398745, 0.80031406, 0]},
              '2018': {'values': [1, 1, 1.00580723, 1.15602798, 1.36748448, 1.55950658, 2.72827293, 1.98275764],
                       'err':    [0, 0, .00179402320, .00563560832, .0232546555, .108132326, .682068233, 1.98275764]}
              }
    

ana_vars = {
    'ak4vars'    : ['Jet_btagDeepB'+LC,'Jet_puId','Jet_jetId'],
    # 'Jet_deepFlavourlepb'+LC, 'Jet_deepFlavouruds'+LC, 'Jet_deepFlavourb'+LC, 'Jet_deepFlavourbb'+LC],
    'ak4lvec'    : {'TLV'         :['JetTLV'+LC],
                    'TLVarsLC'    :['Jet_pt'+LC, 'Jet_eta'+LC, 'Jet_phi'+LC, 'Jet_mass'+LC],
                    'TLVars'      :['Jet_pt', 'Jet_eta', 'Jet_phi', 'Jet_mass'],
                    'jesTotUp'    :['Jet_pt_jesTotalUp', 'Jet_eta' 'Jet_phi', 'Jet_mass_jesTotalUp'],
                    'jesTotDown'  :['Jet_pt_jesTotalDown', 'Jet_eta' 'Jet_phi', 'Jet_mass_jesTotalDown'],
                    'jerUp'       :['Jet_pt_jerUp', 'Jet_eta' 'Jet_phi', 'Jet_mass_jerUp'],
                    'jerDown'     :['Jet_pt_jerDown', 'Jet_eta' 'Jet_phi', 'Jet_mass_jerDown'],
                },
#
    'ak8vars'    : ['FatJet_jetId',
                    'FatJet_deepTagMD_WvsQCD','FatJet_deepTagMD_TvsQCD','FatJet_deepTagMD_bbvsLight',
                    'FatJet_deepTagMD_ZHbbvsQCD',
                    #'FatJet_deepTag_WvsQCD'+LC,'FatJet_deepTag_TvsQCD'+LC,'FatJet_deepTag_ZvsQCD'+LC,
                    #'FatJet_deepTagMD_H4qvsQCD'+LC, 'FatJet_deepTagMD_HbbvsQCD'+LC, 'FatJet_deepTagMD_TvsQCD'+LC, 
                    #'FatJet_deepTagMD_ZbbvsQCD'+LC, 'FatJet_deepTagMD_ZvsQCD'+LC, 'FatJet_deepTagMD_bbvsLight'+LC, 'FatJet_deepTagMD_ccvsLight'+LC,     
                    'FatJet_msoftdrop'+LC,'FatJet_btagDeepB'+LC,'FatJet_btagHbb'+LC,
                    'FatJet_subJetIdx1'+LC,'FatJet_subJetIdx2'+LC],
    'ak8lvec'    : {'TLV'      :['FatJetTLV'+LC],
                    'TLVarsLC' :['FatJet_pt'+LC, 'FatJet_eta'+LC, 'FatJet_phi'+LC, 'FatJet_mass'+LC],
                    'TLVars'   :['FatJet_pt', 'FatJet_eta', 'FatJet_phi', 'FatJet_mass']},
    'ak8sj'      : ['SubJet_pt', 'SubJet_btagDeepB'],
#
    'genpvars'   : ['GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'GenPart_mass', 
                    'GenPart_status', 'GenPart_pdgId', 'GenPart_genPartIdxMother','genTtbarId'], # event level identifier for ttbar+bb
    'genLevCuts' : ['passGenCuts','isZToLL'], # these are MC only
    'valvars'    : ['nJets30'+LC, 'nBottoms'+LC, 
                    #'nSoftBottoms'+LC,
                    #'nResolvedTops'+LC,'nMergedTops'+LC,
                    'passSingleLepElec', 'passSingleLepMu',
                    'MET_phi', 'MET_pt', #'Lep_pt', 'Lep_eta', 'Lep_phi', 'Lep_E',
                    'Pass_IsoTrkVeto', 'Pass_TauVeto', 'Pass_ElecVeto', 'Pass_MuonVeto',
                    'Pass_trigger_muon', 'Pass_trigger_electron'],
    'event'      : ['MET_phi', 'MET_pt','run'],
    'filters'    : ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter',
                    'Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter',
                    'Flag_BadPFMuonFilter','Flag_eeBadScFilter','Flag_ecalBadCalibFilterV2'],
    'HEM_veto'        : ['SAT_Pass_HEMVeto_DataOnly', 'SAT_Pass_HEMVeto_DataAndMC', 'SAT_HEMVetoWeight',
                         'SAT_Pass_HEMVeto_DataOnly'+LC, 'SAT_Pass_HEMVeto_DataAndMC'+LC, 'SAT_HEMVetoWeight'+LC],
    # these are MC only
    'sysvars_mc'      : ['genWeight','puWeight','BTagWeight',
                         'BTagWeight_jes_up','BTagWeight_jes_down',
                         'BTagWeight_lf_up','BTagWeight_lf_down',
                         'BTagWeight_hf_up','BTagWeight_hf_down',
                         'BTagWeight_lfstats1_up','BTagWeight_lfstats1_down',
                         'BTagWeight_lfstats2_up','BTagWeight_lfstats2_down',
                         'BTagWeight_hfstats1_up','BTagWeight_hfstats1_down',
                         'BTagWeight_hfstats2_up','BTagWeight_hfstats2_down',
                         'puWeightUp','puWeightDown', 
                         'pdfWeight_Up','pdfWeight_Down',],
    'sysvars_2016'    : ['PrefireWeight','PrefireWeight_Up','PrefireWeight_Down'],
    'sysvars_2017'    : ['PrefireWeight','PrefireWeight_Up','PrefireWeight_Down'],
    'sysvars_2018'    : [],
    'lheWeights'      : ['PSWeight','LHEScaleWeight','LHEReweightingWeight'],
    'dataHLT_all'     : [ 'HLT_IsoMu24' , 'HLT_IsoMu27', 'HLT_Mu50','HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
                          'HLT_Ele27_WPTight_Gsf', 'HLT_Photon175','HLT_Ele115_CaloIdVT_GsfTrkIdT'],
    'dataHLT_2016'    : ['HLT_IsoTkMu24','HLT_TkMu50','HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50'],
    'dataHLT_2017'    : ['HLT_Ele35_WPTight_Gsf', 'HLT_Ele32_WPTight_Gsf_L1DoubleEG', 'HLT_Photon200', 'HLT_Ele28_eta2p1_WPTight_Gsf_HT150',
                         'HLT_OldMu100','HLT_TkMu100'],
    'dataHLT_2018'    : ['HLT_Ele35_WPTight_Gsf', 'HLT_Ele32_WPTight_Gsf_L1DoubleEG', 'HLT_Photon200', 'HLT_Ele28_eta2p1_WPTight_Gsf_HT150',
                         'HLT_Ele32_WPTight_Gsf',
                         'HLT_OldMu100','HLT_TkMu100'],
    'valRCvars'  : ['ResolvedTopCandidate_discriminator', 'ResolvedTopCandidate_j1Idx', 'ResolvedTopCandidate_j2Idx', 'ResolvedTopCandidate_j3Idx'],
    'label'      : ['isTAllHad']
}

##### DNN backend for Z/H -> bb #####
dnn_ZH_dir  = dataDir+'/NN_files/'
# only event level variables
dnn_ZH_vars = [
    'max_lb_dr','max_lb_invM', 'n_Zh_btag_sj', 'n_ak4jets', 'Zh_score', 'best_rt_score',
    'n_q_outZh', 'n_b_outZh', 'Zh_l_dr', 'n_Zh_sj', 'n_b_inZh', 'Zh_bestb_sj', 'Zh_worstb_sj',
    'Zh_eta','Zh_deepB','b1_outZh_score', 'best_Zh_b_invM_sd', 'Zh_b1_invM_sd', 'Zh_b2_invM_sd','Zh_l_invM_sd',
    'Zh_Wscore', 'Zh_Tscore', 'n_ak8_Zhbb', 'n_ak8jets', 'n_ak4jets',  # repeated?
    'nonZhbb_b1_dr', 'nonZhbb_b2_dr', 
    'Zh_bbscore_sj', 
    'b1_over_Zhpt', 'bb_over_Zhpt',
    'spher','aplan','n_q_inZh']
    #'H_M',
    #'H_pt',
    #'min_lb_invm', 'MET_pt', 'b2oHZpt' 
    #'lb_mtb1', 'lb_invm1', 'lb_dr1',
    ##'weight', 'genWeight'] ####### dont do this...
old_dnn_ZH_vars = [
    'max_lb_dr','max_lb_invm', 'n_H_sj_btag', 'nJets30', 'H_score', 'best_rt_score',
    'n_qnonHbb', 'n_nonHbb', 'Hl_dr', 'n_H_sj', 'n_b_Hbb', 'H_sj_bestb', 'H_sj_worstb',
    'H_eta','H_bbscore','b1_outH_score', 'best_Wb_invM_sd', 'Hb_invM1_sd', 'Hb_invM2_sd','Hl_invm_sd',
    'H_Wscore', 'H_Tscore', 'nhbbFatJets', 'nFatJets', 'nJets', 'nonHbb_b1_dr', 'nonHbb_b2_dr', 
    'H_sj_bbscore', 
    'b1oHZpt', 'bboHZpt',
    'spher','aplan',
    #'H_M',
    #'H_pt',
    #'min_lb_invm', 'MET_pt', 'b2oHZpt' 
    #'lb_mtb1', 'lb_invm1', 'lb_dr1',
    'n_q_Hbb', 'weight', 'genWeight']
uncor_dnn_ZH_vars = [ # tests_failed >= 4 
    'min_lb_dr', 'b2_outH_score'
]
#
dnn_ZH_alpha      = 0.00003 # 200: 0.0003 # 300 : 0.00003
dnn_ZH_batch_size = 512
fl_gamma          = .2 # 200: .1    , 300: 1.5 / .4
fl_alpha          = .85 # 200: .85 , 300: .80 /.85
dnn_ZH_epochs     = 0 #210 ### 200: 120, 300: 100
DNNoutputDir      = dataDir+'/NN_files/'
DNNoutputName     = 'corr_noweight_noM.h5'
DNNmodelName      = 'corr_noweight_model_noM.h5' 
DNNuseWeights     = True
