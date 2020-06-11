import sys
import deepsleepcfg as cfg
from getdata import getData
# testing pbs batch system

getData_cfg = {'files': ['MC_2017'], 'samples': [sys.argv[1]], 'outDir': 'test_ja_files/',
                    'njets':cfg.ZHbbFitCut[1], 'maxJets':cfg.ZHbbFitMaxJets,
                    'treeDir':cfg.tree_dir+'_bb', 'getGenData':True, 'getak8var':True}
#
print('Running getData')
getData(getData_cfg) # in the future this will return ak4, ak8, gen, rtc, and event data
# for now just saves to pkl

