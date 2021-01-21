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
import time
from config.sample_cff import process_cfg, sample_cfg
import config.ana_cff as cfg
#
# add parsargs at some point for year, rundata, minibatch
parser = argparse.ArgumentParser(description='Run analysis over many samples and years using batch system')
parser.add_argument('-s', dest='samples', type=str, 
                    choices=list(process_cfg.keys())+['all','data','mc', 'tt','ttsys'],
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years+['all'],
                    required=True, help='year or all years')
parser.add_argument('-j', dest='jec', type=str, required=False, help='Run with specified jec variation', choices=cfg.jec_variations+['all'], default=None)
parser.add_argument('--script', type=str, required=False, help='Run Analysis or SKim', choices=['runAna','runSkim'], default='runAna')
#parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

log_dir = f'log/{args.year}/'
job_script = f'scripts/{args.script}.sh'

sample_dict = {'all' :process_cfg.keys(),
               'mc'  : [k for k in process_cfg.keys() if 'Data' not in k and 'sys' not in k],
               'data': [k for k in process_cfg.keys() if 'Data' in k],
               'tt'   : process_cfg['TTBar'],
               'ttsys': process_cfg['TTBar_sys'],
               #'tt_bb' : process_cfg['ttbb'],
               #'ttbbsys': process_cfg['ttbb_sys']}
               }

def submit_jobs():
    
    # submit a job for each background / signal / data sample
    if args.samples in sample_dict:
        samples = sample_dict[args.samples]
    else:
        samples = [args.samples]
    years = [args.year] if args.year != 'all' else cfg.Years
    #jecs  = [args.jec]  if args.jec != 'all' else ['JESUp','JESDown','JERUp','JERDown']
    jecs    = cfg.jec_variations if args.jec == 'all' else [args.jec]
    func = submit_runAna if args.script == 'runAna' else submit_runSkim
    #
    for year in years:
        for jec in jecs:
            if jec is not None and re.search(f'201\d', jec) and year not in jec: continue
            if jec and 'HEM' in jec and year != '2018': continue # HEM issue is only for 2018
            func(samples,year,jec)

def submit_runAna(samples, year, jec):
    execute(samples,year,jec)

def submit_runSkim(samples, year, jec):
    for s in samples:
        if 'Data' in s:
            execute_runSkim(process_cfg[s][year],year,jec)
        else:
            execute_runSkim(process_cfg[s],year,jec)

def execute(samples, year, jec):
    for sample in samples:
        add_args  = ''
        add_out_name = ''
        if jec is not None :
            add_args = f',jec={jec}'
            add_out_name = '_'+jec
            tag = jec
        else:
            tag = ''
        #
        if 'Skim' in args.script:
            add_args +=',qsub=True'
        out_name  = sample+'_'+year+add_out_name
        if sample == 'TTBar':
            ppn = 4
        else:
            ppn = 4
        pass_args = f'-v sample={sample},year={year}{add_args}'
        command   = f'qsub -l nodes=1:ppn={ppn} -o {log_dir}Ana_{out_name}std.out -e {log_dir}Ana_{out_name}.stderr '
        command  += f'-N {args.samples}_{args.year}_{sample}{year}{tag} '
        command += f' {pass_args} {job_script}'
        print(command)
        os.system(command)

def execute_runSkim(samples,year,jec):
    popens = []
    for sample in samples:
        args = ['python','runSkim.py','-s',sample,'-y',year,'--qsub','--nopost']
        if jec is not None :
            tag = jec
            args += ['-j',jec]
        else:
            tag = ''
        args += ['-t',tag]
        #popens.append(sb.Popen(args, stdout=sb.PIPE))
        print(args)
        popens.append(sb.Popen(args, stdout=sb.PIPE, stderr=sb.PIPE))
        time.sleep(10)
        num_jobs_running = lambda: int(sb.check_output(
            f"qstat -u $USER -w -f | grep 'Job_Name = runSkim_' | wc -l", shell=True).decode())
        # allow qsub to catch up?
        time.sleep(5)
        while num_jobs_running() > 40:
            time.sleep(30) 
        #popens.append(sb.Popen(args))

if __name__ == '__main__':
    submit_jobs()


