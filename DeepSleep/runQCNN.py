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
                    choices=['nodak8md_dnn_ZH_vars','withbbvl_dnn_ZH_vars','allvars_dnn_ZH_vars'],
                    required=False, help='model input set', default='withbbvl_dnn_ZH_vars')

args = parser.parse_args()

def runQC():
    input_dict = {
        'nodak8md_dnn_ZH_vars': cfg.nodak8md_dnn_ZH_vars,
        'withbbvl_dnn_ZH_vars': cfg.withbbvl_dnn_ZH_vars,
        'allvars_dnn_ZH_vars': cfg.allvars_dnn_ZH_vars,
    }
    for nn_input in input_dict[args.nn_inputs]:
        command = f"qsub -l nodes=1:ppn=4 -N runQCNN_{nn_input} "
        command += f" -o log/nn/{nn_input}.stdout -e log/nn/{nn_input}.stderr "
        command += f'-v nn_input={nn_input} scripts/runQCNN.sh'
        print(command)
        os.system(command)
        #
        num_jobs_running = lambda: int(sb.check_output(
            f"qstat -u $USER -w -f | grep 'Job_Name = runQCNN_' | wc -l", shell=True).decode())
        # allow qsub to catch up?
        time.sleep(5)
        while num_jobs_running() > 30:
            time.sleep(30) 

if __name__ == '__main__':
    runQC()
