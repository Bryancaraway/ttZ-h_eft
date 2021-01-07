import time
import argparse
import re
import config.ana_cff as cfg
from config.sample_cff import sample_cfg
from lib.fun_library import t2Run
#from modules.getdata import getData
from modules.processAna import processAna
from modules.AnaDict    import AnaDict

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('-s', dest='sample', type=str, 
                    #choices=cfg.All_MC+cfg.Data_samples+['test']+cfg.tt_sys_samples+cfg.Sig_EFT_MC+cfg.tt_eft_samples, 
                    choices=sample_cfg.keys(),
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

@t2Run
def runAna ():
    #####
    sample   = args.sample
    isData   = 'Data' in args.sample
    isSignal = re.search(r'TT[Z,H]*\w*', sample) is not None #'TTZH' in sample or 'TTZ_bb' in sample
    isttbar  = sample in cfg.ttbar_samples or sample in cfg.tt_sys_samples or sample in cfg.tt_eft_samples
    tag      = (args.tag + args.jec if args.jec is not None else args.tag)
    #
    input_file = args.inputfile if args.inputfile else cfg.postSkim_dir+f"{args.year}/{sample_cfg[args.sample]['out_name']}/{args.sample}{'.'+tag if tag is not None else ''}.pkl"
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
        ak4_df, ak8_df = [ gD_out[k] for k in ['ak4','ak8'] ]
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


if __name__ == '__main__':
    
    runAna()


