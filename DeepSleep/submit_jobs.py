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
                    choices=list(process_cfg.keys())+['all','data','mc','mc_nosys', 'jec'],
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years+['all'],
                    required=True, help='year or all years')
parser.add_argument('-j', dest='jec', type=str, required=False, help='Run with specified jec variation', choices=cfg.jec_variations+['all'], default=None)
parser.add_argument('--script', type=str, required=False, help='Run Analysis or SKim', choices=['runAna','runSkim','trigSkim'], default='runAna')
parser.add_argument('-q', dest='queue', type=str, required=False, help='Queue to submit jobs to', choices=['hep','econ','james'], default=None)
#parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

log_dir = f'log/{args.year}/'
job_script = f'scripts/{args.script}.sh'

sample_dict = {'all' : [k for k in process_cfg.keys() if 'Data' not in k]+cfg.Data_samples,
               #'mc'  : [k for k in process_cfg.keys() if 'Data' not in k and 'sys' not in k],
               'mc'  : [k for k in process_cfg.keys() if 'Data' not in k],
               'mc_nosys'  : [k for k in process_cfg.keys() if 'Data' not in k and 'sys' not in k],
               'jec' : ['ttZ','ttH','ttbb','single_t','TTBar'],
               'data': cfg.Data_samples,
}

q_limit = {
    'hep'  : 9, 
    'econ' : 20,
    'james': 20,
}

def submit_jobs():
    
    # submit a job for each background / signal / data sample
    if args.samples in sample_dict:
        samples = sample_dict[args.samples]
    else:
        samples = [args.samples]
    years = [args.year] if args.year != 'all' else cfg.Years[::-1]
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
    execute_runAna(samples,year,jec)

def submit_runSkim(samples, year, jec):
    for s in samples:
        if 'Data' in s:
            execute_runSkim(process_cfg[s][year],year,jec)
        else:
            execute_runSkim(process_cfg[s],year,jec)

def execute_runAna(samples, year, jec):
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
        out_name  = sample+'_'+year+add_out_name
        if sample == 'TTBar':
            ppn = 4
        else:
            ppn = 4
        if args.queue is not None:
            add_queue = f'-q {args.queue}'
        else:
            add_queue = ''
        pass_args = f'-v sample={sample},year={year}{add_args}'
        command   = f'qsub {add_queue} -l nodes=1:ppn={ppn} -o {log_dir}Ana_{out_name}.stdout -e {log_dir}Ana_{out_name}.stderr '
        command  += f'-N {args.samples}_{args.year}_{sample}{year}{tag} '
        command += f' {pass_args} {job_script}'
        print(command)
        os.system(command)
        num_jobs_running = lambda: int(sb.check_output(
            #f"qstat -u $USER -w -f | grep 'Job_Name = {args.samples}_{args.year}_{sample}{year}{tag}' | wc -l", shell=True).decode())
            f"qstat -u $USER -w -f | grep 'Job_Name = {args.samples}_{args.year}_{sample}{year}' | wc -l", shell=True).decode())
        # allow qsub to catch up?
        continue_or_wait(num_jobs_running)
        #time.sleep(5)
        #while num_jobs_running() > 40:
        #    time.sleep(30) 

def execute_runSkim(samples,year,jec):
    popens = []
    for sample in samples:
        _args = ['python','runSkim.py','-s',sample,'-y',year,'--qsub','--nopost']
        if args.script == 'trigSkim':
            _args = _args+['--is4trig']
        if jec is not None :
            tag = jec
            _args += ['-j',jec]
        else:
            tag = ''
        if args.queue is not None:
            _args += ['-q',args.queue]
        _args += ['-t',tag]
        #popens.append(sb.Popen(_args, stdout=sb.PIPE))
        print(_args)
        popens.append(sb.Popen(_args, stdout=sb.PIPE, stderr=sb.PIPE))
        time.sleep(10)
        if args.queue is None:
            num_jobs_running = lambda: int(sb.check_output(
                f"qstat -u $USER -w -f | grep 'Job_Name = runSkim_' | wc -l", shell=True).decode())
        else:
            num_jobs_running = (lambda: np.array(sb.check_output(
                f"qstat -q | grep {args.queue}", shell=True).decode().split()[5:7], dtype=int).sum())
        # allow qsub to catch up?
        #print(num_jobs_running())
        continue_or_wait(num_jobs_running)
        #time.sleep(5)
        #while num_jobs_running() > q_limit.get(args.queue,40):
        #    time.sleep(30) 

def continue_or_wait(n_j_running):
    #print(n_j_running())
    time.sleep(5)
    while n_j_running() > q_limit.get(args.queue,60):
        time.sleep(30) 

if __name__ == '__main__':
    submit_jobs()


