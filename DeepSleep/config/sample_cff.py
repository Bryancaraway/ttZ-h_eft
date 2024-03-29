# --- TTZH sample Config File --- #
# python dictionary format

process_cfg = { 
    'Signal_EFT' : ['TTZ_EFT',
                    'TTH_EFT'],
    'Bkg_EFT'    : [#'TTJets_EFT',
                    'TTbb_EFT'],
    
    'ttH'        : ['ttHTobb',
                    'ttHToNonbb'],
    'ttZ'        : ['TTZToLLNuNu',
                    'TTZToBB',    
                    'TTZToQQ'],
    'TTBar'      : ['TTToHadronic',    
                    'TTToSemiLeptonic',
                    'TTTo2L2Nu'],
    'ttbb'       : ['TTbb_Hadronic',
                    'TTbb_SemiLeptonic',
                    'TTbb_2L2Nu'],
    'TTBar_sys'  : ['TTToHadronic_UEDown',   
                    'TTToHadronic_UEUp',     
                    'TTToHadronic_hdampDown',
                    'TTToHadronic_hdampUp',  
                    'TTToSemiLeptonic_UEDown',   
                    'TTToSemiLeptonic_UEUp',     
                    'TTToSemiLeptonic_hdampDown',
                    'TTToSemiLeptonic_hdampUp',  
                    'TTTo2L2Nu_UEDown',   
                    'TTTo2L2Nu_UEUp',     
                    'TTTo2L2Nu_hdampDown',
                    'TTTo2L2Nu_hdampUp'],
    'ttbb_sys'   : ['TTbb_Hadronic_hdampDown',
                    'TTbb_Hadronic_hdampUp',  
                    'TTbb_SemiLeptonic_hdampDown',
                    'TTbb_SemiLeptonic_hdampUp',  
                    'TTbb_2L2Nu_hdampDown',
                    'TTbb_2L2Nu_hdampUp'],
    'single_t'   : ['ST_tW_top',
                    'ST_tW_antitop',
                    'ST_s_lep',
                    'ST_t_top',     
                    'ST_t_antitop',],
    #'tZq_ll',],# might move this to ttX
    #'rare'       : [#'tZq_had',
    #                'THW',
    #                'THQ'],
    'ttX'        : ['TTWJetsToLNu',
                    'TTWJetsToQQ', 
                    'TTTT',  
                    'TTGJets',
                    'THW',
                    'THQ',
                    'tZq_ll',],
    'VJets'      : [#'WJetsToLNu_HT_70to100',
                    #'WJetsToLNu_HT_100to200',  
                    #'WJetsToLNu_HT_200to400',  
                    'WJetsToLNu_HT_400to600',  
                    'WJetsToLNu_HT_600to800',  
                    'WJetsToLNu_HT_800to1200', 
                    'WJetsToLNu_HT_1200to2500',
                    'WJetsToLNu_HT_2500toInf',
                    #'DYJetsToLL_HT_70to100',
                    #'DYJetsToLL_HT_100to200',  
                    #'DYJetsToLL_HT_200to400',  
                    'DYJetsToLL_HT_400to600',  
                    'DYJetsToLL_HT_600to800',  
                    'DYJetsToLL_HT_800to1200',
                    'DYJetsToLL_HT_1200to2500',
                    'DYJetsToLL_HT_2500toInf'],
    'VV'         : ['WW',
                    'WZ',
                    'ZZ'],
    #'VVV'        : ['WWW',
    #                'WWZ',
    #                'WZZ',
    #                'ZZZ',
    #                'WZG',
    #                'WWG'],
    
    'QCD'         : [#'QCD_HT_200to300',
                     #'QCD_HT_300to500',
                     'QCD_HT_500to700',
                     'QCD_HT_700to1000',
                     'QCD_HT_1000to1500',
                     'QCD_HT_1500to2000',
                     'QCD_HT_2000toInf'],
    #
    'Data_SingleElectron' : {'2016':['Data_SingleElectron'],
                             '2017':['Data_SingleElectron'],
                             '2018':['Data_EGamma']},
    'Data_EGamma'         : {'2016':['Data_SingleElectron'],
                             '2017':['Data_SingleElectron'],
                             '2018':['Data_EGamma']},
    'Data_SingleMuon'    :  {'2016':['Data_SingleMuon'],
                             '2017':['Data_SingleMuon'],
                             '2018':['Data_SingleMuon']},
}

signal_xsec = {
    'ttZ': 0.7826,
    'ttH': 0.5084
}

sample_cfg = {                 
    # EFT Samples              
    'TTZ_EFT' : {     'out_name'  : 'TTZ_EFT', 'xs' : 0.7826,  'kf' : 1.0,},                                              
    'TTH_EFT' : {     'out_name'  : 'TTH_EFT', 'xs' : 0.5084,  'kf' : 1.0,},                                              
    'TTJets_EFT' : {  'out_name'  : 'TTJets_EFT',    'xs' : 832.40,  'kf' : 1.0,},                                              
    'TTbb_EFT'   : {  'out_name'  : 'TTbb_EFT',    'xs' : 1.0,     'kf' : 1.0,},                                              
    # ttH                                           
    'ttHTobb'    : {  'out_name'  : 'ttH',        'xs' : .2934,   'kf' : 1.0,},                                              
    'ttHToNonbb' : {  'out_name'  : 'ttH',        'xs' : .215,    'kf' : 1.0,},                                              
    # ttZ                                           
    'TTZToLLNuNu' : { 'out_name'  : 'ttZ',        'xs' : .2529,   'kf' : 1.0,},                                              
    'TTZToBB'     : { 'out_name'  : 'ttZ',        'xs' : .1157,   'kf' : 1.0,},                                              
    'TTZToQQ'     : { 'out_name'  : 'ttZ',        'xs' : .5297,   'kf' : 1.0,},                                              
    # ttbar (5FS)                                   
    'TTToHadronic'     : { 'out_name' : 'TTBar',       'xs' : 380.095, 'kf' : 1.0,},
    'TTToSemiLeptonic' : { 'out_name' : 'TTBar',       'xs' : 364.018, 'kf' : 1.0,},
    'TTTo2L2Nu'        : { 'out_name' : 'TTBar',       'xs' : 88.29,   'kf' : 1.0,},
    # ttbar (4FS)
    'TTbb_Hadronic'     : { 'out_name' : 'ttbb',        'xs' : 1.0,    'kf' : 1.0,},
    'TTbb_SemiLeptonic' : { 'out_name' : 'ttbb',        'xs' : 1.0,    'kf' : 1.0,},
    'TTbb_2L2Nu'        : { 'out_name' : 'ttbb',        'xs' : 1.0,    'kf' : 1.0,},
    # ttbar systematic samples
    'TTToHadronic_UEDown'         : {'out_name' : 'TTBar_sys',     'xs'  : 380.095, 'kf' : 1.0,},                                   
    'TTToHadronic_UEUp'           : {'out_name' : 'TTBar_sys',     'xs'  : 380.095, 'kf' : 1.0,},                                   
    'TTToHadronic_hdampDown'      : {'out_name' : 'TTBar_sys',     'xs'  : 380.095, 'kf': 1.0,},                                   
    'TTToHadronic_hdampUp'        : {'out_name' : 'TTBar_sys',     'xs'  : 380.095, 'kf': 1.0,},
    #
    'TTToSemiLeptonic_UEDown'    : {'out_name' : 'TTBar_sys',      'xs' : 364.018,  'kf' : 1.0,},
    'TTToSemiLeptonic_UEUp'      : {'out_name' : 'TTBar_sys',      'xs' : 364.018,  'kf' : 1.0,},
    'TTToSemiLeptonic_hdampDown' : {'out_name' : 'TTBar_sys',      'xs' : 364.018,  'kf': 1.0,},
    'TTToSemiLeptonic_hdampUp'   : {'out_name' : 'TTBar_sys',      'xs' : 364.018,  'kf': 1.0,},
    #
    'TTTo2L2Nu_UEDown'           : {'out_name' : 'TTBar_sys',      'xs' : 88.29,    'kf' : 1.0,},
    'TTTo2L2Nu_UEUp'             : {'out_name' : 'TTBar_sys',      'xs' : 88.29,    'kf' : 1.0,},
    'TTTo2L2Nu_hdampDown'        : {'out_name' : 'TTBar_sys',      'xs' : 88.29,    'kf': 1.0,},
    'TTTo2L2Nu_hdampUp'          : {'out_name' : 'TTBar_sys',      'xs' : 88.29,    'kf': 1.0,},
    #
    'TTbb_Hadronic_hdampDown'     : {'out_name' : 'ttbb_sys',      'xs' : 1.0,      'kf' : 1.0,},            
    'TTbb_Hadronic_hdampUp'       : {'out_name' : 'ttbb_sys',      'xs' : 1.0,      'kf' : 1.0,},
    #                                                                              
    'TTbb_SemiLeptonic_hdampDown' : {'out_name' : 'ttbb_sys',      'xs' : 1.0,      'kf' : 1.0,},
    'TTbb_SemiLeptonic_hdampUp'   : {'out_name' : 'ttbb_sys',      'xs' : 1.0,      'kf' : 1.0,},
    #                                                                              
    'TTbb_2L2Nu_hdampDown'        : {'out_name' : 'ttbb_sys',      'xs' : 1.0,      'kf' : 1.0,},
    'TTbb_2L2Nu_hdampUp'          : {'out_name' : 'ttbb_sys',      'xs' : 1.0,      'kf' : 1.0,},
    # singleT
    'ST_tW_top'            : {'out_name' : 'single_t',     'xs' : 19.12,   'kf' : 1.0,}, # 35.85 -> 19.12
    'ST_tW_antitop'        : {'out_name' : 'single_t',     'xs' : 19.12,   'kf' : 1.0,}, # 35.85 -> 19.12
    #'ST_tW_top_nofull'     : {'out_name' : 'single_t',     'xs' : 19.12,   'kf' : 1.0,},
    #'ST_tW_antitop_nofull' : {'out_name' : 'single_t',     'xs' : 19.12,   'kf' : 1.0,},
    'ST_s_lep'      : {'out_name' : 'single_t',     'xs' : 6.96,    'kf' : 1.0,},
    'ST_t_top'      : {'out_name' : 'single_t',     'xs' : 136.065, 'kf' : 1.0,},
    'ST_t_antitop'  : {'out_name' : 'single_t',     'xs' : 80.97,   'kf' : 1.0,},
    #'tZq_ll'        : {'out_name' : 'single_t',     'xs' : 0.0758,  'kf' : 1.0,},

    'tZq_had'       : {'out_name' : 'rare',     'xs' : 0.1518,  'kf' : 1.0,}, # not avail for all years
    # ttX
    'tZq_ll'       : {'out_name' : 'ttX',     'xs' : 0.0758,   'kf' : 1.0,},
    'TTWJetsToLNu' : {'out_name' : 'ttX',     'xs' : 0.1793,   'kf' : 1.0,},
    'TTWJetsToQQ'  : {'out_name' : 'ttX',     'xs' : 0.3708,   'kf' : 1.0,},
    'TTTT'         : {'out_name' : 'ttX',     'xs' : 0.009103, 'kf' : 1.0,},
    'TTGJets'      : {'out_name' : 'ttX',     'xs' : 3.697,    'kf' : 1.0,},
    'THW'          : {'out_name' : 'ttX',     'xs' : 0.01517,  'kf' : 1.0,},
    'THQ'          : {'out_name' : 'ttX',     'xs' : 0.07425,  'kf' : 1.0,},
    # Vjets
    'WJetsToLNu_HT_70to100'    : {    'out_name' : 'VJets',    'xs' : 1353,    'kf' : 1.21,},
    'WJetsToLNu_HT_100to200'   : {    'out_name' : 'VJets',    'xs' : 1345.0,  'kf' : 1.21,},
    'WJetsToLNu_HT_200to400'   : {    'out_name' : 'VJets',    'xs' : 359.7,   'kf' : 1.21,},
    'WJetsToLNu_HT_400to600'   : {    'out_name' : 'VJets',    'xs' : 48.91,   'kf' : 1.21,},
    'WJetsToLNu_HT_600to800'   : {    'out_name' : 'VJets',    'xs' : 12.05,   'kf' : 1.21,},
    'WJetsToLNu_HT_800to1200'  : {    'out_name' : 'VJets',    'xs' : 5.501,   'kf' : 1.21,},
    'WJetsToLNu_HT_1200to2500' : {    'out_name' : 'VJets',    'xs' : 1.329,   'kf' : 1.21,},
    'WJetsToLNu_HT_2500toInf'  : {    'out_name' : 'VJets',    'xs' : 0.03216, 'kf' : 1.21,},
    #
    'DYJetsToLL_HT_70to100'     : {   'out_name' : 'VJets',    'xs' : 169.9,    'kf' : 1.23,},
    'DYJetsToLL_HT_100to200'    : {   'out_name' : 'VJets',    'xs' : 147.4,    'kf' : 1.23,},
    'DYJetsToLL_HT_200to400'    : {   'out_name' : 'VJets',    'xs' : 40.99,    'kf' : 1.23,},
    'DYJetsToLL_HT_400to600'    : {   'out_name' : 'VJets',    'xs' : 5.678,    'kf' : 1.23,},
    'DYJetsToLL_HT_600to800'    : {   'out_name' : 'VJets',    'xs' : 1.367,    'kf' : 1.23,},
    'DYJetsToLL_HT_800to1200'   : {   'out_name' : 'VJets',    'xs' : 0.6304,   'kf' : 1.23,},
    'DYJetsToLL_HT_1200to2500'  : {   'out_name' : 'VJets',    'xs' : 0.1514,   'kf' : 1.23,},
    'DYJetsToLL_HT_2500toInf'   : {   'out_name' : 'VJets',    'xs' : 0.003565, 'kf' : 1.23,},
    # VV
    'WW' : {    'out_name' : 'VV',    'xs' :  118.7,     'kf' : 1.0,},
    'WZ' : {    'out_name' : 'VV',    'xs' :  47.13,     'kf' : 1.0,},
    'ZZ' : {    'out_name' : 'VV',    'xs' :  15.8274,   'kf' : 1.0,},
    # VVV
    'WWW' : {   'out_name' : 'VVV',   'xs' :  0.2086,    'kf' : 1.0,},
    'WWZ' : {   'out_name' : 'VVV',   'xs' :  0.1651,    'kf' : 1.0,},
    'WZZ' : {   'out_name' : 'VVV',   'xs' :  0.05565,   'kf' : 1.0,},
    'ZZZ' : {   'out_name' : 'VVV',   'xs' :  0.01398,   'kf' : 1.0,},
    'WZG' : {   'out_name' : 'VVV',   'xs' :  0.04123,   'kf' : 1.0,},
    'WWG' : {   'out_name' : 'VVV',   'xs' :  0.2147,    'kf' : 1.0,},
    # QCD
    #'QCD_HT_200to300'   : {   'out_name' : 'QCD', 'xs' : 1556000, 'kf' : 1.0,},  
    #'QCD_HT_300to500'   : {   'out_name' : 'QCD', 'xs' : 323600,  'kf' : 1.0,},  
    'QCD_HT_500to700'   : {   'out_name' : 'QCD', 'xs' : 29990,   'kf' : 1.0,}, # 'xs2016':32150
    'QCD_HT_700to1000'  : {   'out_name' : 'QCD', 'xs' : 6351,    'kf' : 1.0,}, # 'xs2016':6828
    'QCD_HT_1000to1500' : {   'out_name' : 'QCD', 'xs' : 1093,    'kf' : 1.0,}, # 'xs2016':1208
    'QCD_HT_1500to2000' : {   'out_name' : 'QCD', 'xs' : 99.01,   'kf' : 1.0,}, # 'xs2016':120
    'QCD_HT_2000toInf'  : {   'out_name' : 'QCD', 'xs' : 20.23,   'kf' : 1.0,}, # 'xs2016':25.27
    # Data
    'Data_SingleElectron' : { 'out_name' : 'Data_SingleElectron', 'xs' : 0.0, 'kf' : 1.0,},
    'Data_EGamma'         : { 'out_name' : 'Data_SingleElectron', 'xs' : 0.0, 'kf' : 1.0,},
    'Data_SingleMuon'     : { 'out_name' : 'Data_SingleMuon',     'xs' : 0.0, 'kf' : 1.0,},
}
