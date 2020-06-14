################
## Config for ##
## Deep Sleep ##
################
import operator
#
master_file_path  = './files/'
# Overhead #
file_path         = '/cms/data/store/user/ttxeft/Skim_nanoAOD/'
files             = ['result_2016','result_2017','result_2018_AB','result_2018_CD', 'result_2018']
tree_dir          = 'Training'
MCsamples         = ['TTZ','DY','TTX','DiBoson','TTBarLep']
skim_dir          = master_file_path+'skim/'
skim_kinemFit_dir = master_file_path+'skim_kinemFit/'
skim_Zinv_dir     = master_file_path+'skim_Zinv/'
##
skim_ZHbb_dir     = master_file_path+'skim_ZHbb_200/' # skim_ZHbb
ZHptcut           = 200
##
##############
sample_maxJets  = {'DiLep':{'DY':14, 'TTZ':14, 'TTX':14, 'TTBarLep':11, 'DiBoson':11, 'TriBoson':10},
                   'ZInv':{'WJets':14, 'ZJets':13, 'DiBoson':11, 'TriBoson':11, 'TTX':14, 'QCD':13, 
                           'TTBarHad':13, 'TTBarLep':14, 'TTZ':13 }}
# Kinematic Fit sub cfg args
kinemFitCfg    = (['result_2017'], 
                  [ 'TTX', 'DiBoson', 'TTBarLep', 'TriBoson', 'DY','TTZ'],
                  #['DY'],#['TTZ','DY','TTX', 'DiBoson', 'TTBarLep'], 
                  skim_kinemFit_dir)
kinemFitCut    = (operator.ge, 5)
kinemFitoverlap= 1
kinemFitMaxJets= 14
##### TTZ, Z to MET CONFIG #####
ZinvFitCfg    = (['result_2017'],
                 ['WJets','ZJets','DiBoson','TriBoson','TTX','QCD','TTBarHad','TTBarLep','TTZ'],
                 skim_Zinv_dir)
ZinvFitCut     = (operator.ge, 5)
ZinvFitoverlap = 0
ZinvFitMaxJets = 14
##### TTZ, Z to bb CONFIG #####
ZHbbFitCfg    = (['result_2017'],#
                 #['WJets','ZJets','DY','DiBoson','TriBoson','TTX','QCD','TTBarHad','TTBarLep','TTZ/H'],
                 [ 'TTZH', 'QCD',  'TTX',  'DY', 'WJets', 'TTBarHad',  'DiBoson',  'TriBoson', 'TTBarLep'],#'ZJets'],
                 #[ 'TTZH', 'QCD',  'TTX',  'DY', 'WJets', 'DiBoson',  'TriBoson', 'TTBarLep'],#'ZJets'],
                 
#                 ['TTBarLep'],
                 skim_ZHbb_dir)
ZHbbFitCut    = (operator.ge, 4)
ZHbbFitoverlap = 0
ZHbbFitMaxJets = 100
ZHbb_btagWP    = .4941 # Med for 2017
# ttZ/H->bb SM x-section
ZHbbXsec = {'ttZbb': .1157,
            'ttHbb': .2934 }
ZHbbtotXsec = ZHbbXsec['ttZbb'] + ZHbbXsec['ttHbb']
# ttZ/H->bb MC count 2017
n_ZHbbMC_dict      = {'ttZbb': 163876,
                      'ttHbb': 5698653 }
n_ZHbbMC           = n_ZHbbMC_dict['ttZbb'] + n_ZHbbMC_dict['ttHbb']
# Train Overhead #
train_dir      = master_file_path+'train/'
train_over_dir = master_file_path+'train_overSample/'
test_dir       = master_file_path+'test/'
val_dir        = master_file_path+'val/'

###################
# Input Variables #
LC = '_drLeptonCleaned'
#
ana_vars = {
    'ak4vars'    : ['Jet_btagDeepB'+LC],
    'ak4lvec'    : {'TLV'         :['JetTLV'+LC],
                    'TLVarsLC'    :['Jet_pt'+LC, 'Jet_eta'+LC, 'Jet_phi'+LC, 'Jet_mass'+LC],
                    'TLVars'      :['Jet_pt', 'Jet_eta', 'Jet_phi', 'Jet_mass'],
                    'jesTotUp'    :['Jet_pt_jesTotalUp', 'Jet_eta' 'Jet_phi', 'Jet_mass_jesTotalUp'],
                    'jesTotDown'  :['Jet_pt_jesTotalDown', 'Jet_eta' 'Jet_phi', 'Jet_mass_jesTotalDown'],
                    'jerUp'       :['Jet_pt_jerUp', 'Jet_eta' 'Jet_phi', 'Jet_mass_jerUp'],
                    'jerDown'     :['Jet_pt_jerDown', 'Jet_eta' 'Jet_phi', 'Jet_mass_jerDown'],
                },
#
    'ak8vars'    : ['FatJet_tau1'+LC,'FatJet_tau2'+LC,'FatJet_tau3'+LC,'FatJet_tau4'+LC,
                    'FatJet_deepTag_WvsQCD'+LC,'FatJet_deepTag_TvsQCD'+LC,'FatJet_deepTag_ZvsQCD'+LC,
                    'FatJet_msoftdrop'+LC,'FatJet_btagDeepB'+LC,'FatJet_btagHbb'+LC,
                    'FatJet_subJetIdx1'+LC,'FatJet_subJetIdx2'+LC],
    'ak8lvec'    : {'TLV'      :['FatJetTLV'+LC],
                    'TLVarsLC' :['FatJet_pt'+LC, 'FatJet_eta'+LC, 'FatJet_phi'+LC, 'FatJet_mass'+LC],
                    'TLVars'   :['FatJet_pt', 'FatJet_eta', 'FatJet_phi', 'FatJet_mass']},
    'ak8sj'      : ['SubJet_pt', 'SubJet_btagDeepB'],
#
    'genpvars'   : ['GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'GenPart_E', 'GenPart_status', 'GenPart_pdgId', 'GenPart_genPartIdxMother'], # these are MC only
    'genLevCuts' : ['passGenCuts','isZToLL'], # these are MC only
    'valvars'    : ['nResolvedTops'+LC,'nMergedTops'+LC,'nBottoms'+LC,'nSoftBottoms'+LC,'nJets30'+LC,
                    'MET_phi', 'MET_pt', 'Lep_pt', 'Lep_eta', 'Lep_phi', 'Lep_E',
                    'Pass_IsoTrkVeto', 'Pass_TauVeto', 'Pass_ElecVeto', 'Pass_MuonVeto',
                    'Stop0l_trigger_eff_Electron_pt', 'Stop0l_trigger_eff_Muon_pt', 
                    'Pass_trigger_muon', 'Pass_trigger_electron'],
    'sysvars'    : ['genWeight','weight','BTagWeight','puWeight','ISRWeight','PrefireWeight', # these are MC only
                    'BTagWeight_Up', 'BTagWeight_Down', 'puWeight_Up','puWeight_Down', 'pdfWeight_Up','pdfWeight_Down',
                   'ISRWeight_Up','ISRWeight_Down','PrefireWeight_Up','PrefireWeight_Down'],
    'valRCvars'  : ['ResolvedTopCandidate_discriminator', 'ResolvedTopCandidate_j1Idx', 'ResolvedTopCandidate_j2Idx', 'ResolvedTopCandidate_j3Idx'],
    'label'      : ['isTAllHad']
}

# Derived Varialbes #
ak4comb = 'true'
ak8comb = 'true'
# Model Hyperparams #
epochs        = 300
alpha         = 0.001
batch_size    = int(32768/4)
hiddenl       = 3 
NNoutputDir   = './NN_ouput/'
NNoutputName  = 'sixth_try.h5' 
NNmodel1Name  = 'sixth_model.h5'
#NNOUTPUTNAME  = 'fifth_try.h5' 
#NNmodel1Name  = 'fifth_model.h5'
useWeights    = True
plotHist      = True

##### CNN backend prep #####
skim_cnn_dir  = master_file_path+'skim_cnn/'
cnn_data_dir  = (master_file_path+'train_cnn/', master_file_path+'test_cnn/', master_file_path+'val_cnn/')
cnnCut        = (operator.ge, 5)
cnnMaxJets    = 10
cnnProcessCfg = (files, MCsamples, skim_cnn_dir)
cnn_vars = ['btagCSVV2','btagDeepB', 'qgl', 'pt', 'eta', 'phi', 'E']
#
cnn_alpha      = 0.001
cnn_batch_size = 512
cnn_epochs     = 30
CNNoutputDir   = './CNN_output/'
CNNoutputName  = 'first_try.h5'
CNNmodelName   = 'first_model.h5'

##### DNN backend for Z/H -> bb #####
dnn_ZH_dir  = (skim_ZHbb_dir+'train_ZH/', skim_ZHbb_dir+'test_ZH/', skim_ZHbb_dir+'val_ZH/')
aux_ZH_dir  = skim_ZHbb_dir+'aux_ZH/'
# only event level variables
dnn_ZH_vars = [
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
DNNoutputDir      = skim_ZHbb_dir+'DNN_ZH_output/'
DNNoutputName     = 'corr_noweight_noM.h5'
DNNmodelName      = 'corr_noweight_model_noM.h5' 
DNNuseWeights     = True
