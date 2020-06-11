# run full analysis on kodiak's batch system
# use case: python submitjobs.py
#
import sys
import os
import deepsleepcfg as cfg
#

# submit a job for each background / signal / data sample
samples = cfg.ZHbbFitCfg[1]
log_dir = 'log/'
for sample in samples:
    command = f'qsub -o {log_dir}{sample}.out -e {log_dir}{sample}.err -v sample={sample} test.sh'
    print(command)
    os.system(command)

