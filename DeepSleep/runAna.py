import sys
import time
import argparse
import cfg.deepsleepcfg as cfg
from modules.getdata import getData
from modules.processAna import processAna

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('-s', dest='sample', type=str, choices=cfg.MC_samples+cfg.Data_samples+cfg.Pow_samples, 
                    required=True, help='sample to analyze')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years,
                    required=True, help='year')
parser.add_argument('-i', dest='roofile', type=str, required=False, help="Optional input root file, leave out '.root'", default=None)
parser.add_argument('-t', dest='tag', type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

class runAna ():
    start = time.perf_counter()
    #####
    sample   = args.sample # should change to parsargs  at some point
    roofile  = (args.roofile if args.roofile is not None else (f'Data_{args.year}' if 'Data' in args.sample else f'MC_{args.year}') ) 
    isData   = 'Data' in roofile
    isSignal = 'TTZH'     in sample
    isttbar  = 'TTBarLep' in sample

    #####

    print('Running getData...')
    getData_cfg = {'roofile': roofile, 'sample': sample, 'outDir': 'files/', 'year':args.year,
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
    processAna_cfg = {'outDir': 'files/', 'outTag':args.tag, 'year':args.year, 'isData':isData, 'isSignal':isSignal, 'isttbar':isttbar,
                      'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df,
                      'sample':sample}
    processAna(processAna_cfg)

    #####
    finish = time.perf_counter()
    print(f'\nTime to finish {__name__} for {args.sample}: {finish-start:.1f}\n')

if __name__ == '__main__':
    
    runAna()


