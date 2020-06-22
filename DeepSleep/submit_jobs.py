# run full analysis on kodiak's batch system
# use case: python submitjobs.py
#
import sys
import os
import argparse
import cfg.deepsleepcfg as cfg
#
# add parsargs at some point for year, rundata, minibatch
parser = argparse.ArgumentParser(description='Run analysis over many samples and years using batch system')
parser.add_argument('-s', dest='samples', type=str, choices=['all','mc','data'], 
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years,
                    required=True, help='year')
args = parser.parse_args()

sample_dict = {'all':cfg.MC_samples+cfg.Data_samples,
               'mc':cfg.MC_samples,
               'data':cfg.Data_samples}
#cfg.MC_samples+cfg.Data_samples
# submit a job for each background / signal / data sample
job_script = 'scripts/runAna.sh'
samples = sample_dict[args.samples]
log_dir = 'log/'
os.system('rm log/*')
for sample in samples:
    command = f'qsub -l nodes=1:ppn=1 -o {log_dir}{sample}_{args.year}.out -e {log_dir}{sample}_{args.year}.err -v sample={sample},year={args.year} {job_script}'
    print(command)
    os.system(command)

