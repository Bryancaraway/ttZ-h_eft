import time
import argparse
import re
import config.ana_cff as cfg
from config.sample_cff import sample_cfg, process_cfg
from lib.fun_library import t2Run
#from modules.getdata import getData
from modules.processAna import processAna
from modules.AnaDict    import AnaDict

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('-s', dest='sample', type=str, 
                    #choices=cfg.All_MC+cfg.Data_samples+['test']+cfg.tt_sys_samples+cfg.Sig_EFT_MC+cfg.tt_eft_samples, 
                    choices=list(sample_cfg.keys())+list(process_cfg.keys()),
                    required=True, help='sample to analyze')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years,
                    required=True, help='year')
parser.add_argument('-i', dest='inputfile', type=str, required=False, help="Optional input pkl file", default=None)
parser.add_argument('--jec', dest='jec',     type=str, required=False, help='Run with specified jec variation', choices=cfg.jec_variations+[''], default=None)
parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
parser.add_argument('--estart', dest='estart', type=int, required=False, help='parse event to start from', default=None)
parser.add_argument('--estop',  dest='estop',  type=int, required=False, help='parse event to stop at', default=None)
parser.add_argument('--condor', action='store_true', required=False, help='Flag is running on Condor', default=False)
parser.add_argument('--keep_all', action='store_true', required=False, help='Flag to keep ak4,ak8,gen,rtc arrays', default=False)
args = parser.parse_args()

if args.jec is not None and re.search(r'201\d', str(args.jec)) and args.year not in args.jec:
    raise ValueError(f"{args.jec} is not compatible with year choice: {args.year}")
    exit()

def runAna():
    #if 'Data' in args.sample:
    #    for sample in process_cfg[args.sample][args.year]:
    #        analyzer(sample)
    if  args.sample in sample_cfg.keys() or 'Data' in args.sample:
        analyzer(args.sample)
    elif args.sample in process_cfg.keys():
        for sample in process_cfg[args.sample]:
            print(sample)
            analyzer(sample)
        # run post job
        post_job(process_cfg[args.sample])

@t2Run
def analyzer(sample):
    #####
    sample   = sample
    isData   = 'Data' in sample
    #isSignal = re.search(r'TT[Z,H]*\w*', sample) is not None #'TTZH' in sample or 'TTZ_bb' in sample
    isSignal = re.search(r'(tt(Z|H))', sample_cfg[sample]['out_name']) is not None
    isttbar  = re.search(r'(TT|tt)(Bar)|(bb)', sample_cfg[sample]['out_name']) is not None
    #isttbar  = re.search(r'TT[To,bb]', sample) is not None
    print(isSignal, isttbar)
    #isttbar  = sample in cfg.ttbar_samples or sample in cfg.tt_sys_samples or sample in cfg.tt_eft_samples
    tag      = (args.tag + args.jec if args.jec is not None else args.tag)
    #
    input_file = args.inputfile if args.inputfile else cfg.postSkim_dir+f"{args.year}/{sample_cfg[sample]['out_name']}/{sample}{'.'+tag if tag != '' else ''}.pkl"
    if isData and (args.jec is not None and args.jec != ''): exit()
    #####

    #print('Running getData...')
    #getData_cfg = {'roofile': roofile, 'sample': sample, 'outDir': 'files/', 'year':args.year,
    #               'njets':cfg.ZHbbFitMinJets, 'maxAk4Jets':cfg.ZHbbFitMaxJets,
    #               'treeDir':cfg.tree_dir+'_bb', 'isData':isData, 'jec_sys': args.jec, 'jec_type': args.jetjec,
    #               'estart':args.estart, 'estop':args.estop}
    #gD_out = getData(getData_cfg).getdata()
    gD_out = AnaDict.read_pickle(f'{input_file}')
    if isData: 
        ak4_df, ak8_df = [ AnaDict(gD_out[k]) for k in ['ak4','ak8'] ]
        val_df = gD_out['events']
        gen_df = None
    else:
        ak4_df, ak8_df, gen_df = [ AnaDict(gD_out[k]) for k in ['ak4','ak8','gen'] ]
        val_df = gD_out['events']

    #####
    print('Running processAna...')
    processAna_cfg = {'outDir': 'files/', 'outTag':tag, 'year':args.year, 'isData':isData, 'isSignal':isSignal, 'isttbar':isttbar,
                      'ak4_df':ak4_df, 'ak8_df':ak8_df , 'val_df':val_df, 'gen_df':gen_df,
                      'sample':sample, 'condor':args.condor, 'keep_all': args.keep_all}
    processAna(processAna_cfg)

    #####

def post_job(samples):
    import pandas as pd
    import os
    out_df = pd.DataFrame()
    #out_name = sample_cfg[samples[0]]['out_name']
    wDir   = f"{cfg.master_file_path}/{args.year}/{'mc_files' if 'Data' not in args.sample else 'data_files'}/"
    outTag = '_'+args.jec if args.jec is not None else ''
    for sample in samples:
        target_sample = f'{wDir}/{sample}{outTag}_val.pkl'
        target_df = pd.read_pickle(target_sample)
        out_df = pd.concat([out_df, target_df], axis='rows', ignore_index=True)
        os.system(f'rm {target_sample}')
    out_df.to_pickle(f'{wDir}/{args.sample}{outTag}_val.pkl')

if __name__ == '__main__':
    
    runAna()


