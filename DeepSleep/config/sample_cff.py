# --- TTZH sample Config File --- #
# python dictionary format


sample_cfg = {
    # EFT Samples
    'TTZ_EFT' : {
        'out_name'  : 'Signal_EFT',
        'xs'        : 0.7826,
        'kf'        : 1.0,
    },
    'TTH_EFT' : {
        'out_name'  : 'Signal_EFT',
        'xs'        : 0.5084,
        'kf'        : 1.0,
    },
    'TTJets_EFT' : {
        'out_name'  : 'Bkg_EFT',
        'xs'        : 832.40,
        'kf'        : 1.0,
    },
    'TTbb_EFT' : {
        'out_name'  : 'Bkg_EFT',
        'xs'        : 1.0,
        'kf'        : 1.0,
    },
    # ttH
    'ttHTobb' : {
        'out_name'  : 'ttH',
        'xs'        : .2934,
        'kf'        : 1.0,
    },
    'ttHToNonbb' : {
        'out_name'  : 'ttH',
        'xs'        : .215,
        'kf'        : 1.0,
    },
    # ttZ
    'TTZToLLNuNu' : {
        'out_name'  : 'ttZ',
        'xs'        : .2529,
        'kf'        : 1.0,
    },
    'TTZToBB' : {
        'out_name'  : 'ttZ',
        'xs'        : .1157,
        'kf'        : 1.0,
    },
    'TTZToQQ' : {
        'out_name'  : 'ttZ',
        'xs'        : .5297,
        'kf'        : 1.0,
    },
    # ttbar (5FS)
    'TTToHadronic' : {
        'out_name' : 'TTBar',
        'xs'       : 380.095,
        'kf'       : 1.0,
    },
    'TTToSemiLeptonic' : {
        'out_name' : 'TTBar',
        'xs'       : 364.017888,
        'kf'       : 1.0,
    },
    'TTTo2L2Nu' : {
        'out_name' : 'TTBar',
        'xs'       : 88.29,
        'kf'       : 1.0,
    },
    # ttbar (4FS)
    'TTbb_Hadronic' : {
        'out_name' : 'ttbb',
        'xs'       : 1.0,
        'kf'       : 1.0,
    },
    'TTbb_SemiLeptonic' : {
        'out_name' : 'ttbb',
        'xs'       : 1.0,
        'kf'       : 1.0,
    },
    'TTbb_2L2Nu' : {
        'out_name' : 'ttbb',
        'xs'       : 1.0,
        'kf'       : 1.0,
    },
    # ttbar systematic samples
    # ttbar (5FS)
    # ttbar (4FS)
    # singleT
    # ttX
    # Vjets
    # VV
    # VVV
    # Data
}
