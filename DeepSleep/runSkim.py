import os
import time
import argparse
import re
import numpy as np
from glob import glob
#
import config.ana_cff as cfg
from config.sample_cff import sample_cfg
from lib.fun_library import fillne, t2Run
#  module classes
from modules.AnaDict import AnaDict
from modules.preSkim import PreSkim
from modules.Skim    import Skim

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('-s', dest='sample', type=str, 
                    choices=sample_cfg.keys(),
                    #cfg.All_MC+cfg.Data_samples+['test']+cfg.tt_sys_samples+cfg.Sig_EFT_MC+cfg.tt_eft_samples, 
                    required=True, help='sample to analyze')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years,
                    required=True, help='year')
#parser.add_argument('-i', dest='roofile', type=str, required=False, help="Optional input root file, leave out '.root'", default=None)
parser.add_argument('--jec', dest='jec',     type=str, required=False, help='Run with specified jec variation', choices=cfg.jec_variations+[''], default=None)
parser.add_argument('--jjec', dest='jetjec', type=str, required=False, help='Which jet to compute jec', choices=['ak4','ak8','ak8jmsr'], default=None)
parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
parser.add_argument('--qsub', dest='qsub',  action='store_true', required=False, help='Run jobs on pbs', default=False)
#parser.add_argument('--estart', dest='estart', type=int, required=False, help='parse event to start from', default=None)
#parser.add_argument('--estop',  dest='estop',  type=int, required=False, help='parse event to stop at', default=None)
#parser.add_argument('--condor', action='store_true', required=False, help='Flag is running on Condor', default=False)
#parser.add_argument('--keep_all', action='store_true', required=False, help='Flag to keep ak4,ak8,gen,rtc arrays', default=False)
args = parser.parse_args()


@t2Run
def runSkim():
    sample_dir = cfg.postproc_dir
    ####
    isData = 'Data' in args.sample
    if not isData: # is MC
        files = glob(f"{sample_dir}/{args.year}/{args.sample}_{args.year}/*.root")
    else : # is data
        files = glob(f"{sample_dir}/{args.year}/{args.sample}_{args.year}_Period*/*.root")
    isttbar = 'TTTo' in args.sample or 'TTJets' in args.sample
    isttbb  = 'TTbb' in args.sample
    out_dir = f"{cfg.postSkim_dir}/{args.year}/{sample_cfg[args.sample]['out_name']}"
    if not os.path.exists(out_dir):
        os.system(f'mkdir {out_dir}')
    #### Pre Skim
    print('Running PreSkim')
    run_preSkim = PreSkim(files, args.sample, args.year, isttbar=isttbar, isttbb=isttbb)
    metaData    = run_preSkim.get_metadata()
    #### Skim
    print('Running Skim')
    for i, sfile in enumerate(files):
        print(sfile)
        out_file = f'{out_dir}/{args.sample}_{i}.pkl'
        run_Skim = Skim(sfile, args.sample, args.year, isData)
        Skim_data = run_Skim.get_skim()
        #print(Skim_data)
        #Skim_data['metaData'] = metaData
        AnaDict(Skim_data).to_pickle(out_file)
    #### Post Skim
    pkl_files = glob(f'{out_dir}/{args.sample}_*.pkl')
    final_pkl = {}
    for i, pkl in enumerate(pkl_files):
        pkl_dict = AnaDict.read_pickle(pkl)
        if i == 0:
            final_pkl = pkl_dict
            os.system(f'rm {pkl}')
            continue
        print(len(final_pkl['ak8']['FatJet_pt']))
        for k in pkl_dict:
            #print(pkl_dict[k])
            for var in pkl_dict[k]:
                try:
                    final_pkl[k][var] = concatenate([final_pkl[k][var],pkl_dict[k][var]])
                except AttributeError:
                    final_pkl[k][var] = np.concatenate([final_pkl[k][var],pkl_dict[k][var]])
            #
        os.system(f'rm {pkl}')
        #
    #
    final_pkl['metaData'] = metaData
    AnaDict(final_pkl).to_pickle(f'{out_dir}/{args.sample}.pkl')
        
def concatenate(jaggedarrays):
    import awkward
    contents = np.concatenate([j.flatten() for j in jaggedarrays])
    counts = np.concatenate([j.counts for j in jaggedarrays])
    return awkward.JaggedArray.fromcounts(counts, contents)

if __name__ == '__main__':
    runSkim()
