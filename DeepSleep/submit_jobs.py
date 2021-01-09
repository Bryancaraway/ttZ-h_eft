# run full analysis on kodiak's batch system
# use case: python submitjobs.py
#
import sys
import os
import argparse
import numpy as np
import pandas as pd
import re
import subprocess as sb
from config.sample_cff import process_cfg, sample_cfg
import config.ana_cff as cfg
#
# add parsargs at some point for year, rundata, minibatch
parser = argparse.ArgumentParser(description='Run analysis over many samples and years using batch system')
parser.add_argument('-s', dest='samples', type=str, 
                    choices=list(process_cfg.keys())+['all','data','mc', 'tt','ttsys','ttbb','ttbbsys'],
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years+['all'],
                    required=True, help='year or all years')
parser.add_argument('-j', dest='jec', type=str, required=False, help='Run with specified jec variation', choices=cfg.jec_variations+['all'], default=None)
parser.add_argument('--script', type=str, required=False, help='Run Analysis or SKim', choices=['runAna','runSkim'], default='runAna')
#parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

def submit_jobs():

    sample_dict = {'all' :process_cfg.keys(),
                   'mc'  : [k for k in process_cfg.keys() if 'Data' not in k and 'sys' not in k],
                   'data': [k for k in process_cfg.keys() if 'Data' in k],
                   'tt'   : process_cfg['TTBar'],
                   'ttsys': process_cfg['TTBar_sys'],
                   'ttbb' : process_cfg['ttbb'],
                   'ttbbsys': process_cfg['ttbb_sys']}
    
    # submit a job for each background / signal / data sample
    if args.samples in process_cfg:
        samples = process_cfg[args.samples]
    else:
        samples = sample_dict.get(args.samples,[args.samples])
    years = [args.year] if args.year != 'all' else cfg.Years
    #jecs  = [args.jec]  if args.jec != 'all' else ['JESUp','JESDown','JERUp','JERDown']
    jecs    = cfg.jec_variations if args.jec == 'all' else [args.jec]
    #
    iterate_args(execute, samples, years, jecs)
    if 'Skim' in args.script: # no need to go further
        return
    #
    num_jobs_running = lambda: int(sb.check_output(
        f"qstat -u $USER -w -a | grep {args.samples}_{args.year}_ | wc -l", shell=True).decode())
    while num_jobs_running() > 0:
        time.sleep(30)
    #
    iterate_args(postjob, samples, years, jecs)


def iterate_args(func, samples, years, jecs):
    for year in years:
        #os.system(f'rm {log_dir}*{year}*')
        for jec in jecs:
            if jec is not None and re.search(f'201\d', jec) and year not in jec:
                continue
            func(samples, year, jec)
                

def execute(samples, year, jec):
    log_dir = f'log/{year}/'
    job_script = f'scripts/{args.script}.sh'
    for sample in samples:
        add_args  = ''
        add_out_name = ''
        if jec is not None:
            add_args = f',jec={jec}'
            add_out_name = '_'+jec
            #
        if 'Skim' in args.script:
            add_args +=',qsub=True'
        out_name  = sample+'_'+year+add_out_name
        if sample == 'TTBar' and year == '2018':
            ppn = 4
        else:
            ppn = 1
        pass_args = f'-v sample={sample},year={year}{add_args}'
        command   = f'qsub -l nodes=1:ppn={ppn} -o {log_dir}{out_name}.out -e {log_dir}{out_name}.err '
        command  += f'-N {args.samples}_{args.year}_{sample}{year}{jec} {pass_args} {job_script}'
        print(command)
        os.system(command)

def postjob(samples, year, jec, out_df=None, wDir=None):
    out_df = pd.DataFrame(), 
    wDir   = f"{cfg.master_file_path}/{year}/{'mc_files' if 'Data' not in args.samples else 'data_files'}/"
    outTag = '_'+jec if jec is not None else ''
    for sample in samples:
        target_sample = f'{wDir}/{sample}{outTag}_val.pkl'
        if len(out_df) == 0:
            out_df = pd.read_pickle(target_sample)
        else:
            out_df = pd.concat([out_df,pd.read_pickle(target_sample)], axis='rows', ignore_axis=True)
    out_df.to_pickle(f'{wDir}/{args.samples}{outTag}_val.pkl')
    

if __name__ == '__main__':
    submit_jobs()

