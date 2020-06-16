# run full analysis on kodiak's batch system
# use case: python submitjobs.py
#
import sys
import os
import cfg.deepsleepcfg as cfg
#
# add parsargs at some point for year, rundata, minibatch
# submit a job for each background / signal / data sample
job_script = 'scripts/runAna.sh'
samples = cfg.ZHbbFitCfg[1]
log_dir = 'log/'
os.system('rm log/*')
for sample in samples:
    command = f'qsub -l nodes=1:ppn=1 -o {log_dir}{sample}.out -e {log_dir}{sample}.err -v sample={sample} {job_script}'
    print(command)
    os.system(command)

