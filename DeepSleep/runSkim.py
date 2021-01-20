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
from config.sample_cff import sample_cfg, process_cfg
from lib.fun_library import t2Run
#  module classes
from modules.AnaDict  import AnaDict
from modules.Skim     import Skim
from modules.PostSkim import PostSkim

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
parser.add_argument('--noskim', dest='noskim',  action='store_true', required=False, help='Run postSkim', default=False)
parser.add_argument('--nopost', dest='nopost',  action='store_true', required=False, help='Run postSkim', default=False)
args = parser.parse_args()

if args.jec is not None and re.search(r'201\d', str(args.jec)) and args.year not in args.jec:
    raise ValueError(f"{args.jec} is not compatible with year choice: {args.year}")
#
isData = 'Data' in args.sample
sample_dir = cfg.postproc_dir
log_dir = f'log/{args.year}/{args.sample}/'
if not os.path.exists(log_dir):
    os.system(f'mkdir {log_dir}')
#
tag     = ('.'+(f'{args.tag}_' if args.tag != '' else '')+f'{args.jec}' if args.jec is not None else args.tag)

@t2Run
def runSkim():
    #
    files = get_jobfiles()
    #
    print(f"{sample_dir}/{args.year}/{args.sample}_{args.year}/*.root")
    isttbar = 'TTTo' in args.sample or 'TTJets' in args.sample
    isttbb  = 'TTbb' in args.sample
    out_dir = f"{cfg.postSkim_dir}/{args.year}/{sample_cfg[args.sample]['out_name']}"
    print(out_dir)
    if not os.path.exists(out_dir):
        os.system(f'mkdir {out_dir}')
    #### Skim #######
    print(args.nopost)
    print('Running Skim')
    if not args.noskim:
        if args.qsub:
            parallel_skim(files, out_dir, tag)
        else:
            sequential_skim(files, out_dir, tag)
    #### Post Skim ####
    if not args.nopost:
        postSkim = PostSkim(args.sample, args.year, isData, out_dir, tag)
        postSkim.run()
        #postSkim(out_dir, tag)
        
def sequential_skim(files, out_dir, tag):
    golden_json=json.load(open(cfg.goodLumis_file[args.year]))
    import multiprocessing
    from functools import partial
    pool = multiprocessing.Pool(8)
    #for i, sfile in enumerate(files):
    _worker = partial( skim_worker, out_dir=out_dir,tag=tag,golden_json=golden_json)  
    _ = pool.map(_worker, zip(range(0,len(files)),files))
    pool.close()

def skim_worker(i_sfile,out_dir,tag,golden_json):
        print(i_sfile)
        i, sfile = i_sfile
        out_file = f'{out_dir}/'+(args.outfile.replace('.pkl',f'_{i}.pkl') if args.outfile else f'{args.sample}_{i}{tag}.pkl')
        #
        run_Skim = Skim(sfile, args.sample, args.year, isData, jec_sys=args.jec,  golden_json=golden_json)
        Skim_data = run_Skim.get_skim()
        #
        AnaDict(Skim_data).to_pickle(out_file)
        del run_Skim, Skim_data

def parallel_skim(files, out_dir, tag):
    # use glob
    os.system(f'rm {out_dir}/*{args.sample}*{tag}*.pkl ; rm {log_dir}Skim_{args.sample}_*{tag}*.std*') # start of job, get rid of old files
    job_script = 'scripts/runSkim.sh'
    for i, sfile in enumerate(files):
        # rerun runSkim without pre/post skim using pbs
        command = f"qsub -l nodes=1:ppn=8 -N runSkim_{args.sample}_{args.year}{tag}_{i}  "
        command += f" -o {log_dir}Skim_{args.sample}_{tag}{i}.stdout -e {log_dir}Skim_{args.sample}_{tag}{i}.stderr "
        add_args  = ''
        if args.jec is not None:
            add_args = f',jec={args.jec}'
            #
        out_name  = f'{args.sample}_{i}{tag}.pkl'
        pass_args = f'-v sample={args.sample},year={args.year},infile={sfile},outfile={out_name},nopost=True{add_args}'
        command += f'{pass_args} {job_script}'
        print(command)
        os.system(command)
    # make sure jobs are finished before exiting
    num_jobs_running = lambda: int(sb.check_output(
            f"qstat -u $USER -w -f | grep 'Job_Name = runSkim_{args.sample}_{args.year}{tag}_' | wc -l", shell=True).decode())
    # allow qsub to catch up?
    time.sleep(5)
    print(f"qstat -u $USER -w -f | grep 'Job_Name = runSkim_{args.sample}_{args.year}{tag}_' | wc -l", num_jobs_running())
    while num_jobs_running() > 0:
        time.sleep(30) 
    # jobs are finished here
    # run postjob
    command = f"qsub -l nodes=1:ppn=1 -N PostSkim_{args.sample}_{args.year}{tag} "
    command += f" -o {log_dir}{args.sample}{tag}.stdout -e {log_dir}{args.sample}{tag}.stderr "
    if args.jec is not None:
        add_args = f',jec={args.jec}'
        #
    pass_args = f'-v sample={args.sample},year={args.year},noskim=True{add_args}'
    command += f'{pass_args} {job_script}'
    print(command)
    os.system(command)
    

def get_jobfiles():
    files = []
    #if not isData: # is MC
    #    pre_skim_files = glob(f"{sample_dir}/{args.year}/{args.sample}_{args.year}/*.root")
    #else : # is data
    #    pre_skim_files = glob(f"{sample_dir}/{args.year}/{args.sample}_{args.year}_Period*/*.root")
    #
    if args.roofile:
        if '.list' in args.roofile:
            with open(args.roofile, 'r') as lf:
                files = [l.strip('\n') for l in lf.readlines()]
            os.system(f'rm {args.roofile}') # clear up directory
        else:
            files = [args.roofile]
    else:
        file_header = args.sample if not isData else process_cfg[args.sample][args.year][0]
        files = glob(f"{sample_dir}/{args.year}/{file_header}_{args.year}"+("_Period*" if isData else "")+"/*.root")
    #
    if len(files) > 10 and args.qsub and '.list' not in files[0]: # create filelist of size 5 for less strain on kodiak
        new_files = []
        size = 48 if len(files) > 300 else 12
        os.system(f'rm {log_dir}{args.sample}_*{tag}*.list') # get rid of old list files
        for i in range(0,len(files),size):
            new_file = log_dir+f'{args.sample}_{tag}{i//size}.list'
            new_files.append(new_file)
            with open(new_file, 'w') as lf:
                try:
                    lf.writelines([f+'\n' for f in files[i:i+size]])
                except:
                    lf.writelines([f+'\n' for f in files[i:]])
        files = new_files
    return files

if __name__ == '__main__':
    runSkim()
