import sys
import time
import cfg.deepsleepcfg as cfg
from modules.getdata import getData
from modules.processAna import processAna

class runAna ():
    start = time.perf_counter()
    #####
    sample   = sys.argv[1] # should change to parsargs  at some point
    roofile  = ('Data_2017' if 'Data' in sys.argv[1] else 'MC_2017')
    isData   = 'Data' in roofile
    isSignal = 'TTZH'     in sample
    isttbar  = 'TTBarLep' in sample
    #
    print('Running getData...')
    getData_cfg = {'roofile': roofile, 'sample': sys.argv[1], 'outDir': 'files/',
                   'njets':cfg.ZHbbFitMinJets, 'maxAk4Jets':cfg.ZHbbFitMaxJets,
                   'treeDir':cfg.tree_dir+'_bb', 'isData':isData}

    gD_out = getData(getData_cfg).getdata()
    if isData: 
        ak4_df, ak8_df, val_df, rtc_df = gD_out
        gen_df = None
    else:
        ak4_df, ak8_df, val_df, rtc_df, gen_df = gD_out
        #####
    print('Running processAna...')
    processAna_cfg = {'outDir': 'files/', 'year':'2017', 'isData':isData, 'isSignal':isSignal, 'isttbar':isttbar,
                      'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df,
                      'sample':sys.argv[1]}
    processAna(processAna_cfg)
    #####
    finish = time.perf_counter()
    print(f'\nTime to finish {__name__} for {sys.argv[1]}: {finish-start:.1f}\n')

if __name__ == '__main__':
    
    runAna()


