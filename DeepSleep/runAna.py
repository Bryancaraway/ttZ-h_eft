import time
import argparse
import config.ana_cff as cfg
from modules.getdata import getData
from modules.processAna import processAna

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('-s', dest='sample', type=str, choices=cfg.All_MC+cfg.Data_samples+['test']+cfg.tt_sys_samples, 
                    required=True, help='sample to analyze')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years,
                    required=True, help='year')
parser.add_argument('-i', dest='roofile', type=str, required=False, help="Optional input root file, leave out '.root'", default=None)
parser.add_argument('--jec', dest='jec',     type=str, required=False, help='Run with specified jec variation', choices=['JESUp','JESDown','JERUp','JERDown',''], default=None)
parser.add_argument('--jjec', dest='jetjec', type=str, required=False, help='Which jet to compute jec', choices=['ak4','ak8'], default=None)
parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
parser.add_argument('--estart', dest='estart', type=int, required=False, help='parse event to start from', default=None)
parser.add_argument('--estop',  dest='estop',  type=int, required=False, help='parse event to stop at', default=None)
parser.add_argument('--condor', action='store_true', required=False, help='Flag is running on Condor', default=False)
parser.add_argument('--keep_all', action='store_true', required=False, help='Flag to keep ak4,ak8,gen,rtc arrays', default=False)
args = parser.parse_args()

class runAna ():
    start = time.perf_counter()
    #####
    sample   = args.sample if args.sample != 'test' else 'TriBoson' # should change to parsargs  at some point
    roofile  = (args.roofile if args.roofile is not None else (f'Data_{args.year}' if 'Data' in args.sample else f'MC_{args.year}') ) 
    isData   = 'Data' in roofile
    isSignal = 'TTZH' in sample or 'TTZ_bb' in sample
    print(isSignal)
    isttbar  = sample in cfg.ttbar_samples or sample in cfg.tt_sys_samples
    if (args.jetjec is None and args.jec is not None) or (args.jetjec is not None and args.jec is None):
        raise('Must use --jjec AND --jec to run analysis with jec variation ')
    tag      = (args.tag + args.jetjec +args.jec if args.jec is not None else args.tag)
    #tag      = (tag + agrs.estop if args.estop is not None)
    if isData and (args.jec is not None and args.jec != ''): exit()
    #####

    print('Running getData...')
    getData_cfg = {'roofile': roofile, 'sample': sample, 'outDir': 'files/', 'year':args.year,
                   'njets':cfg.ZHbbFitMinJets, 'maxAk4Jets':cfg.ZHbbFitMaxJets,
                   'treeDir':cfg.tree_dir+'_bb', 'isData':isData, 'jec_sys': args.jec, 'jec_type': args.jetjec,
                   'estart':args.estart, 'estop':args.estop}
    gD_out = getData(getData_cfg).getdata()
    if isData: 
        ak4_df, ak8_df, val_df, rtc_df = gD_out
        gen_df = None
    else:
        ak4_df, ak8_df, val_df, rtc_df, gen_df = gD_out

    #####
    print('Running processAna...')
    processAna_cfg = {'outDir': 'files/', 'outTag':tag, 'year':args.year, 'isData':isData, 'isSignal':isSignal, 'isttbar':isttbar,
                      'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df, 'rtc_df':rtc_df,
                      'sample':sample, 'condor':args.condor, 'keep_all': args.keep_all}
    processAna(processAna_cfg)

    #####
    finish = time.perf_counter()
    print(f'\nTime to finish {__name__} for {args.sample}: {finish-start:.1f}\n')

if __name__ == '__main__':
    
    runAna()


