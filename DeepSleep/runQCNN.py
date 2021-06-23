# job sub script for qc datacards

import os
import time
import argparse
import re
import numpy as np
import subprocess as sb
#
import config.ana_cff as cfg
#  module classes
from modules.AnaDict     import AnaDict


parser = argparse.ArgumentParser(description='Run qc datacard script for nn inputs')
parser.add_argument('--nn_inputs', dest='nn_inputs', type=str, 
                    choices=['nodak8md_dnn_ZH_vars','withbbvl_dnn_ZH_vars','allvars_dnn_ZH_vars','hl2'],
                    required=False, help='model input set', default='withbbvl_dnn_ZHgenm_vars')
parser.add_argument('--rerun', action='store_true', required=False, help='Rerun files which failed (expert use only)', default=False)

args = parser.parse_args()

dc_dir = 'Higgs-Combine-Tool/'

def runQC():
    input_dict = {
        #'nodak8md_dnn_ZH_vars': cfg.nodak8md_dnn_ZH_vars,
        'withbbvl_dnn_ZH_vars': cfg.withbbvl_dnn_ZH_vars,
        'withbbvl_dnn_ZHgenm_vars' : cfg.withbbvl_dnn_ZHgenm_vars,
        'hl2': [f'NN_{i}' for i in np.arange(0,64)],
        #'allvars_dnn_ZH_vars': cfg.allvars_dnn_ZH_vars,
    }
    for nn_input in input_dict[args.nn_inputs]:
        if args.rerun and check_for_file(nn_input) == True:
            continue
        command = f"qsub -l nodes=1:ppn=8 -N runQCNN_{nn_input} "
        #command = f"qsub -q hep -l nodes=gpu006:ppn=8 -N runQCNN_{nn_input} "
        command += f" -o log/nn/{nn_input}.stdout -e log/nn/{nn_input}.stderr "
        command += f'-v nn_input={nn_input} scripts/runQCNN.sh'
        print(command)
        os.system(command)
        #
        num_jobs_running = lambda: int(sb.check_output(
            f"qstat -u $USER -w -f | grep 'Job_Name = runQCNN_' | wc -l", shell=True).decode())
        # allow qsub to catch up?
        time.sleep(1)
        while num_jobs_running() > 60:
            time.sleep(30) 

def check_for_file(nn_input):
    return any([True for f in os.listdir(dc_dir) if nn_input in f])


if __name__ == '__main__':
    runQC()
