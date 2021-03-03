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
#### NanoAODv7 PostProcessed Sample Directory ####
preproc_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PreProcessed/'
postproc_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/'
postSkim_dir = '/cms/data/store/user/ttxeft/NanoAODv7/Skim/'
# Overhead #
import os
if   os.path.exists('/cms/data/store/user/ttxeft/') : # test to see if on kodiak
    file_path         = '/cms/data/store/user/ttxeft/Skim_nanoAOD/' # for kodiak
elif os.path.exists('/eos/uscms/') or 'condor' in _cdir: # test to see if on lpc will need to fix for condor on kodiak i think
    file_path        = 'root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/'
    preproc_dir  = preproc_dir.replace('/cms/data','/eos/uscms').replace('ttxeft','bcaraway')
    postproc_dir = postproc_dir.replace('/cms/data','/eos/uscms').replace('ttxeft','bcaraway')
    postSkim_dir = postSkim_dir.replace('/cms/data','/eos/uscms').replace('ttxeft','bcaraway')
else: raise("Not on Kodiak or LPC, please manually input file_path in file: ./config/ana_cff.py")

##
nn                = 'withbbvl_NN'
ZHptcut           = 200
Years             = ['2016','2017','2018']
## EFT ##
Sig_EFT_MC        = ['TTZ_EFT','TTH_EFT']
tt_eft_samples    = ['TTJets_EFT','TTBB_EFT']
#########
Sig_MC            = ['ttH','ttZ']
Bkg_MC            = ['TTBar','ttbb','ttX','single_t','VV','VVV','VJets']
All_MC            = ['ttZ','ttH','TTBar','ttbb','single_t','ttX','VV','VVV','VJets']

# Handle systematic sample docs
#tt_sys_samples    = ['TTBar_UEUp','TTBar_UEDown','TTBar_hdampUp','TTBar_hdampDown',
#                     'ttbb_hdampUp','ttbb_hdampDown']
tt_sys_samples    = ['TTBar_sys','ttbb_sys']
tt_bb             = ['ttbb',]
tt_bb_sys         = [ 'ttbb_hdampUp','ttbb_hdampDown', ]

#
#jec_variations    = [jtype+jec for jec in ['JESUp','JESDown','JERUp','JERDown'] for jtype in ['ak4','ak8']]
jecs              = [jec+y for jec in ['jesRelativeSample','jesHF' , 'jesAbsolute', 'jesEC2', 'jesBBEC1'] for y in Years] + \
                    ['jesHF' , 'jesAbsolute', 'jesEC2', 'jesBBEC1', 'jesRelativeBal', 'jesFlavorQCD',
                     'jesHEMIssue',  # 2018 only 
                     'ak4jer','ak8jer', # jer
                     'jms','jmr'] # puppi sdm corr
jec_variations    = [jec+ud for jec in jecs for ud in ['Up','Down']]

sig_sys_samples   = [sig+'_'+jec for sig in Sig_MC for jec in jec_variations]
jec_bkg           = ['TTBar','ttbb','single_t'] # for the sake of memory
bkg_sys_samples   = [bkg+'_'+jec for bkg in jec_bkg for jec in jec_variations] + tt_sys_samples
all_sys_samples   = sig_sys_samples + bkg_sys_samples
#
Data_samples      = ['Data_SingleElectron','Data_SingleMuon']

Lumi              = {'2016': 35.917149,
                     '2017': 41.525338,
                     '2018': 59.72444,
                     '2018preHEM' : 21.1,
                     '2018postHEM': 38.6,
                     'run2': 137.166648,
                     'Total': 137.166648
                  } 
goodLumis_file   = {
    '2016':dataDir+'/good_lumis/'+'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt', # re-reco
    '2017':dataDir+'/good_lumis/'+'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt', # re-reco
    '2018':dataDir+'/good_lumis/'+'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt', # re-reco
    #'2016':dataDir+'/good_lumis/'+'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',  # legacy
    #'2017':dataDir+'/good_lumis/'+'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',# legacy
    #'2018':dataDir+'/good_lumis/'+'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',  # legacy
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
                                    #(x['Muon_FlagId'] >= 1)    &  # not in NanoAODv7 files
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

ana_vars = {
    'ak4vars'    : ['Jet_btagDeepB','Jet_puId','Jet_jetId',],
    'ak4mcvars'  : ['Jet_btagSF_deepcsv_shape',
                    'Jet_btagSF_deepcsv_shape_up_jes','Jet_btagSF_deepcsv_shape_down_jes',
                    'Jet_btagSF_deepcsv_shape_up_hf','Jet_btagSF_deepcsv_shape_down_hf',
                    'Jet_btagSF_deepcsv_shape_up_lf','Jet_btagSF_deepcsv_shape_down_lf',
                    'Jet_btagSF_deepcsv_shape_up_hfstats1','Jet_btagSF_deepcsv_shape_down_hfstats1',
                    'Jet_btagSF_deepcsv_shape_up_lfstats1','Jet_btagSF_deepcsv_shape_down_lfstats1',
                    'Jet_btagSF_deepcsv_shape_up_hfstats2','Jet_btagSF_deepcsv_shape_down_hfstats2',
                    'Jet_btagSF_deepcsv_shape_up_lfstats2','Jet_btagSF_deepcsv_shape_down_lfstats2',
                    'Jet_btagSF_deepcsv_shape_up_cferr1','Jet_btagSF_deepcsv_shape_down_cferr1',
                    'Jet_btagSF_deepcsv_shape_up_cferr2','Jet_btagSF_deepcsv_shape_down_cferr2',
                ],
    # 'Jet_deepFlavourlepb'+LC, 'Jet_deepFlavouruds'+LC, 'Jet_deepFlavourb'+LC, 'Jet_deepFlavourbb'+LC],
    'ak4lvec'    : {'TLVarsLC'    :['Jet_pt'+LC, 'Jet_eta'+LC, 'Jet_phi'+LC, 'Jet_mass'+LC],
                    'TLVars'      :['Jet_pt', 'Jet_eta', 'Jet_phi', 'Jet_mass'],
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
    'ak8lvec'    : {'TLVarsLC' :['FatJet_pt'+LC, 'FatJet_eta'+LC, 'FatJet_phi'+LC, 'FatJet_mass'+LC],
                    'TLVars'   :['FatJet_pt', 'FatJet_eta', 'FatJet_phi', 'FatJet_mass']},
    'ak8sj'      : ['SubJet_pt', 'SubJet_btagDeepB'],
#
    'genpvars'   : ['GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'GenPart_mass', 
                    'GenPart_status', 'GenPart_pdgId', 'GenPart_genPartIdxMother','genTtbarId'], # event level identifier for ttbar+bb

    'event'      : ['MET_phi', 'MET_pt',
                    'MET_T1_phi', 'MET_T1_pt',
                    'MET_T1Smear_phi', 'MET_T1Smear_pt',
                    'PV_npvsGood','run','luminosityBlock','event'],
    'filters_all'    : ['Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter',
                        'Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter',
                        'Flag_BadPFMuonFilter','Flag_eeBadScFilter'],
    'filters_year' : {'2016': [], '2017':['Flag_ecalBadCalibFilterV2'], '2018':['Flag_ecalBadCalibFilterV2']},
    # these are MC only
    'sysvars_mc'      : ['genWeight','puWeight',
                         #'BTagWeight',
                         #'BTagWeight_jes_up','BTagWeight_jes_down',
                         #'BTagWeight_lf_up','BTagWeight_lf_down',
                         #'BTagWeight_hf_up','BTagWeight_hf_down',
                         #'BTagWeight_lfstats1_up','BTagWeight_lfstats1_down',
                         #'BTagWeight_lfstats2_up','BTagWeight_lfstats2_down',
                         #'BTagWeight_hfstats1_up','BTagWeight_hfstats1_down',
                         #'BTagWeight_hfstats2_up','BTagWeight_hfstats2_down',
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
nodak8md_dnn_ZH_vars = [
    'outZH_b1_pt','outZH_b2_pt',
    'outZH_b1_score','outZH_b2_score',
    'outZh_q1_pt','outZh_q2_pt',
    'outZh_q1_btag','outZh_q2_btag',
    #
    'outZH_b1_q_mindr','outZH_b2_q_mindr',
    'outZH_q_q_dr_nearb1','outZH_q_q_dr_nearb2',
    'outZH_qq_M_nearb1','outZH_qq_M_nearb2',
    'outZH_b1q_M',#'outZH_b2q_M',
    'outZH_b1_qq_dr','outZH_b2_qq_dr',
    'outZH_b1qq_M','outZH_b2qq_M',
    'ZH_b1qq_dr','ZH_b2qq_dr',
    'ZH_lbb1qq_dr','ZH_lbb2qq_dr',
    #'l_b1_mtb',
    'l_b2_mtb',
    #
    #'HT',
    'Zh_closeb_invM',#'Zh_closeq_invM',
    'n_ak8jets', 'n_ak4jets','n_ak8_Zhbb',
    #'Zh_Wscore', 'Zh_Tscore', 
    #'outZh_max_Wscore', 'outZh_max_Tscore',
    #'Zh_bbvLscore', 'outZh_max_bbvLscore',
    #'outZh_max_ak8pt',
    'outZh_max_ak8sdM',
    'outZh_b12_m', 'outZh_b12_dr', 
    'ht_b', 'ht_outZh',
    #
    'n_Zh_btag_sj',
    'Zh_bestb_sj', 
    'Zh_worstb_sj',
    #'sjpt1_over_Zhpt','sjpt2_over_Zhpt',
    #
    'nonZhbb_q1_dr', 
    'nonZhbb_b1_dr',
    'inZhb_outZhb_dr',
    #
    'Zh_l_dr', 'Zh_l_invM_sd', 
    'l_b1_invM','l_b2_invM',
    'l_b1_dr','l_b2_dr',
    #
    'spher','aplan',
    'n_b_inZh', 'n_q_inZh',
    'n_b_outZh', 'n_q_outZh'
]
withbbvl_dnn_ZH_vars = nodak8md_dnn_ZH_vars+['Zh_bbvLscore']

allvars_dnn_ZH_vars = [
    'outZH_b1_pt','outZH_b2_pt',
    'outZH_b1_score','outZH_b2_score',
    'outZh_q1_pt','outZh_q2_pt',
    'outZh_q1_btag','outZh_q2_btag',
    #
    'outZH_b1_q_mindr','outZH_b2_q_mindr',
    'outZH_q_q_dr_nearb1','outZH_q_q_dr_nearb2',
    'outZH_qq_M_nearb1','outZH_qq_M_nearb2',
    'outZH_b1q_M',#'outZH_b2q_M',
    'outZH_b1_qq_dr','outZH_b2_qq_dr',
    'outZH_b1qq_M','outZH_b2qq_M',
    'ZH_b1qq_dr','ZH_b2qq_dr',
    'ZH_lbb1qq_dr','ZH_lbb2qq_dr',
    #'l_b1_mtb',
    'l_b2_mtb',
    #
    #'HT',
    'Zh_closeb_invM',#'Zh_closeq_invM',
    'n_ak8jets', 'n_ak4jets','n_ak8_Zhbb',
    'Zh_Wscore', 'Zh_Tscore', 
    'outZh_max_Wscore', 'outZh_max_Tscore',
    'Zh_bbvLscore', 'outZh_max_bbvLscore',
    #'outZh_max_ak8pt',
    'outZh_max_ak8sdM',
    'outZh_b12_m', 'outZh_b12_dr', 
    'ht_b', 'ht_outZh',
    #
    'n_Zh_btag_sj',
    'Zh_bestb_sj', 
    'Zh_worstb_sj',
    #'sjpt1_over_Zhpt','sjpt2_over_Zhpt',
    #
    'nonZhbb_q1_dr', 
    'nonZhbb_b1_dr',
    'inZhb_outZhb_dr',
    #
    'Zh_l_dr', 'Zh_l_invM_sd', 
    'l_b1_invM','l_b2_invM',
    'l_b1_dr','l_b2_dr',
    #
    'spher','aplan',
    'n_b_inZh', 'n_q_inZh',
    'n_b_outZh', 'n_q_outZh'
]

selvars_dnn_ZH_vars = [
    'outZH_b1_pt','outZH_b2_pt',
    'outZH_b1_score',#'outZH_b2_score',
    #'outZh_q1_pt','outZh_q2_pt',
    'outZh_q1_btag',#'outZh_q2_btag',
    #
    'outZH_b1_q_mindr','outZH_b2_q_mindr',
    'outZH_q_q_dr_nearb1',#'outZH_q_q_dr_nearb2',
    #'outZH_qq_M_nearb1','outZH_qq_M_nearb2',
    'outZH_b1q_M',#'outZH_b2q_M',
    #'outZH_b1_qq_dr','outZH_b2_qq_dr',
    #'outZH_b1qq_M','outZH_b2qq_M',
    'ZH_b1qq_dr',#'ZH_b2qq_dr',
    'ZH_lbb1qq_dr','ZH_lbb2qq_dr',
    #'l_b1_mtb',
    'l_b2_mtb',
    #
    #'HT',
    'Zh_closeb_invM',#'Zh_closeq_invM',
    'n_ak8jets', 'n_ak4jets','n_ak8_Zhbb',
    'Zh_Wscore', #'Zh_Tscore', 
    'outZh_max_Wscore', 'outZh_max_Tscore',
    'Zh_bbvLscore', 'outZh_max_bbvLscore',
    #'outZh_max_ak8pt',
    'outZh_max_ak8sdM',
    'outZh_b12_m', 'outZh_b12_dr', 
    'ht_b', #'ht_outZh',
    #
    'n_Zh_btag_sj',
    #'Zh_bestb_sj', 
    'Zh_worstb_sj',
    #'sjpt1_over_Zhpt','sjpt2_over_Zhpt',
    #
    #'nonZhbb_q1_dr', 
    'nonZhbb_b1_dr',
    'inZhb_outZhb_dr',
    #
    'Zh_l_dr', 'Zh_l_invM_sd', 
    'l_b1_invM','l_b2_invM',
    'l_b1_dr','l_b2_dr',
    #
    'spher','aplan',
    'n_b_inZh', 'n_q_inZh',
    'n_b_outZh',# 'n_q_outZh',
]

old_dnn_ZH_vars = [
    'max_lb_dr',
    #'max_lb_invM','best_Zh_b_invM_sd',
    'n_Zh_btag_sj', 'Zh_bbvLscore', 'outZh_max_bbvLscore',#'best_rt_score',
    'n_b_outZh','n_Zh_sj',  'Zh_bestb_sj', #'Zh_worstb_sj',
    #'Zh_eta', 'Zh_l_dr','n_q_outZh', 
    'Zh_deepB','b1_outZh_score', 'b2_outZh_score', 'Zh_b1_invM_sd', 'Zh_b2_invM_sd','Zh_l_invM_sd',
    'Zh_Wscore', 'Zh_Tscore', 'outZh_max_Wscore', 'outZh_max_Tscore', 
    'max_farl_b_q_dr',
    'ht_b', 'ht_outZh', 'min_farl_b_q_dr', 'outZh_bqq_mass', 
    'outZh_bb_dr', 'outZh_qq_dr',
    #'Zh_bqq_dr', 
    'Zh_lbbqq_dr',
    #'n_ak8_Zhbb', 
    'n_ak8jets', 'n_ak4jets',  
    'nonZhbb_b1_dr', 'nonZhbb_b2_dr', 
    #'sjpt1_over_Zhpt', 'sjpt2_over_Zhpt',
    #'Zh_bbscore_sj', 
    #'b1_over_Zhpt', 'bb_over_Zhpt',
    'spher','aplan','n_b_inZh', 'n_q_inZh']
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
