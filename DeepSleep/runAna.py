import sys
import time
import deepsleepcfg as cfg
from getdata import getData
from processAna import processAna
# testing pbs batch system
start = time.perf_counter()
#####
sample   = sys.argv[1] # should change to parsargs  at some point
roofile  = 'MC_2017'
isData   = 'Data' in roofile
isSignal = 'TTZH'     in sample
isttbar  = 'TTBarLep' in sample

getData_cfg = {'roofile': roofile, 'samples': [sys.argv[1]], 'outDir': 'files/',
               'njets':cfg.ZHbbFitCut[1], 'maxJets':cfg.ZHbbFitMaxJets,
               'treeDir':cfg.tree_dir+'_bb', 'getGenData':True, 'getak8var':True}
#
print('Running getData...')
ak4_df, ak8_df, val_df, gen_df, rtc_df = getData(getData_cfg).getdata()

print('Running processAna...')
processAna_cfg = {'outDir': 'files/', 'year':'2017', 'isData':isData, 'isSignal':isSignal, 'isttbar':isttbar,
                  'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df,
                  'sample':sys.argv[1]}
processAna(processAna_cfg)
#####
finish = time.perf_counter()
print(f'\nTime to finish {__name__} for {sys.argv[1]}: {finish-start:.1f}\n')


