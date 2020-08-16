# run full analysis on kodiak's batch system
# use case: python submitjobs.py
#
import sys
import os
import argparse
import config.ana_cff as cfg
#
# add parsargs at some point for year, rundata, minibatch
parser = argparse.ArgumentParser(description='Run analysis over many samples and years using batch system')
parser.add_argument('-s', dest='samples', type=str, choices=['all','mc','data'], 
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years+['all'],
                    required=True, help='year or all years')
#parser.add_argument('-j', dest='jec',     type=str, required=False, help='Run with specified jec variation', choices=['JESUp','JESDown','JERUp','JERDown','all'], default='')
#parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

sample_dict = {'all':cfg.All_MC+cfg.Data_samples,
               'mc':cfg.All_MC,
               'data':cfg.Data_samples}
#cfg.MC_samples+cfg.Data_samples
# submit a job for each background / signal / data sample
job_script = 'scripts/runAna.sh'
samples = sample_dict[args.samples]
years = [args.year] if args.year != 'all' else cfg.Years
#jecs  = [args.jec]  if args.jec != 'all' else ['JESUp','JESDown','JERUp','JERDown']

for year in years:
    log_dir = f'log/'
    os.system(f'rm log/*{year}*')
    for sample in samples:
        #for jec in jecs:
         #   if jec == '':
        #    command = f'qsub -l nodes=1:ppn=1 -o {log_dir}{sample}{jec}_{year}.out -e {log_dir}{sample}{jec}_{year}.err -v sample={sample},year={year} {job_script}'
        #else:
        command = f'qsub -l nodes=1:ppn=1 -o {log_dir}{sample}_{year}.out -e {log_dir}{sample}_{year}.err -v sample={sample},year={year} {job_script}'
        print(command)
        os.system(command)

