import os
import time
import argparse
import re
import numpy as np
import pandas as pd
import json
from glob import glob
import subprocess as sb
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
parser.add_argument('-i', dest='roofile', type=str, required=False, help="Optional input root file", default=None)
parser.add_argument('-o', dest='outfile', type=str, required=False, help="Optional output name", default=None)
parser.add_argument('-j','--jec', dest='jec', type=str, required=False, help='Run with specified jec variation', choices=cfg.jec_variations+[''], default=None)
parser.add_argument('-t', dest='tag', type=str, required=False, help='Optional tag to add to output file', default='')
parser.add_argument('--qsub', dest='qsub',  action='store_true', required=False, help='Run jobs on pbs', default=False)
parser.add_argument('--nopre', dest='nopre',  action='store_true', required=False, help='Run preSkim', default=False)
parser.add_argument('--nopost', dest='nopost',  action='store_true', required=False, help='Run postSkim', default=False)
#parser.add_argument('--estart', dest='estart', type=int, required=False, help='parse event to start from', default=None)
#parser.add_argument('--estop',  dest='estop',  type=int, required=False, help='parse event to stop at', default=None)
#parser.add_argument('--condor', action='store_true', required=False, help='Flag is running on Condor', default=False)
#parser.add_argument('--keep_all', action='store_true', required=False, help='Flag to keep ak4,ak8,gen,rtc arrays', default=False)
args = parser.parse_args()

isData = 'Data' in args.sample
if args.jec is not None and re.search(r'201\d', str(args.jec)) and args.year not in args.jec:
        raise ValueError(f"{args.jec} is not compatible with year choice: {args.year}")
        exit()

@t2Run
def runSkim():

    sample_dir = cfg.postproc_dir
    ####
    if args.roofile:
        files = [args.roofile]
    else:
        if not isData: # is MC
            files = glob(f"{sample_dir}/{args.year}/{args.sample}_{args.year}/*.root")
        else : # is data
            files = glob(f"{sample_dir}/{args.year}/{args.sample}_{args.year}_Period*/*.root")
    #
    print(f"{sample_dir}/{args.year}/{args.sample}_{args.year}/*.root")
    isttbar = 'TTTo' in args.sample or 'TTJets' in args.sample
    isttbb  = 'TTbb' in args.sample
    out_dir = f"{cfg.postSkim_dir}/{args.year}/{sample_cfg[args.sample]['out_name']}"
    tag     = ('.'+(f'{args.tag}_' if args.tag != '' else '')+f'{args.jec}' if args.jec is not None else args.tag)
    if not os.path.exists(out_dir):
        os.system(f'mkdir {out_dir}')
    #### Pre Skim ######
    print('Running PreSkim')
    metaData = {}
    if not args.nopre:
        #
        run_preSkim = PreSkim(files, args.sample, args.year, isttbar=isttbar, isttbb=isttbb)
        metaData    = run_preSkim.get_metadata()
        #
    #### Skim #######
    print('Running Skim')
    if args.qsub:
        parallel_skim(files, out_dir, tag)
    else:
        sequential_skim(files, out_dir, tag)
    #### Post Skim ####
    if not args.nopost:
        postSkim(metaData, out_dir, tag)
        
def sequential_skim(files, out_dir, tag):
    golden_json=json.load(open(cfg.goodLumis_file[args.year]))
    for i, sfile in enumerate(files):
        print(sfile)
        out_file = f'{out_dir}/'+(args.outfile if args.outfile else f'{args.sample}_{i}{tag}.pkl')
        #
        run_Skim = Skim(sfile, args.sample, args.year, isData, jec_sys=args.jec,  golden_json=golden_json)
        Skim_data = run_Skim.get_skim()
        #
        AnaDict(Skim_data).to_pickle(out_file)
        del run_Skim, Skim_data

def parallel_skim(files, out_dir, tag):
    log_dir = f'log/{args.year}/'
    job_script = 'scripts/runSkim.sh'
    for i, sfile in enumerate(files):
        # rerun runSkim without pre/post skim using pbs
        command = f"qsub -l nodes=1:ppn=1 -n {args.sample}_{args.year}{tag}_{i}  "
        command += " -o {log_dir}{args.sample} -e {log_dir}{args.sample} "
        add_args  = ''
        if args.jec is not None:
            add_args = f',jec={args.jec}'
            #
        out_name  = f'{args.sample}_{i}{tag}.pkl'
        pass_args = f'-v sample={args.sample},year={args.year},infile={sfile},outfile={out_name}.pkl,noprepost=True,{add_args}'
        command += f'{pass_args} {job_script}'
        print(command)
        os.system(command)
    # make sure jobs are finished before exiting
    num_jobs_running = lambda: int(sb.check_output(
            "qstat -u $USER | grep {args.sample}_{args.year}{tag} | wc -l", shell=True).decode())
    while num_jobs_running() > 0:
        time.sleep(30)
    # jobs are finished here


def postSkim(metaData, out_dir, tag):
    # todo : add weight, update btag weight with metadata, sample name to events
    pkl_files = glob(f'{out_dir}/{args.sample}_*{tag}.pkl')
    final_pkl = {}
    for i, pkl in enumerate(pkl_files):
        pkl_dict = AnaDict.read_pickle(pkl)
        if i == 0:
            final_pkl = pkl_dict
            os.system(f'rm {pkl}')
            continue
        #print(len(final_pkl['ak8']['FatJet_pt']))
        for k in pkl_dict:
            #print(pkl_dict[k])
            if k == 'events' :
                print(len(final_pkl[k]))
                final_pkl[k] = pd.concat([final_pkl[k],pkl_dict[k]], axis='rows', ignore_index=True)
                print(len(final_pkl[k]))
                continue
            for var in pkl_dict[k]:
                try:
                    final_pkl[k][var] = concatenate([final_pkl[k][var],pkl_dict[k][var]]) # jagged
                except AttributeError:
                    final_pkl[k][var] = np.concatenate([final_pkl[k][var],pkl_dict[k][var]]) # might not work with pandas dataframe
            #
        os.system(f'rm {pkl}')
        #
    #
    final_pkl['metaData'] = metaData
    #add normalization, sample name, apply r factor to BtagWeights
    final_pkl['events']['weight'] = (metaData['xs']*metaData['kf']*cfg.Lumi[args.year]*1000)/metaData['tot_events']
    final_pkl['events']['sample'] = metaData['sample']
    # to do btagweight
    final_pkl['events'].loc[:,'BTagWeight'] = final_pkl['events']['BTagWeight']*metaData['r_nom']
    for btagw in final_pkl['events'].filter(like='BTagWeight_', axis='columns').keys():
        final_pkl['events'].loc[:,btagw] = final_pkl['events'][btagw]*metaData['r_'+btagw.replace('BTagWeight_','')]
    #
    AnaDict(final_pkl).to_pickle(f'{out_dir}/{args.sample}{tag}.pkl')

def concatenate(jaggedarrays):
    import awkward
    contents = np.concatenate([j.flatten() for j in jaggedarrays])
    counts = np.concatenate([j.counts for j in jaggedarrays])
    return awkward.JaggedArray.fromcounts(counts, contents)

if __name__ == '__main__':
    runSkim()
