'''
Script used to transfer files from temporary 
condor job directory to main working directory 
file system
'''

import sys
import subprocess as sb
import os
import re
import argparse
if __name__ == '__main__':
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
    # this should mean that wherever the script is executed, it should work
    os.system(f'cd {sys.path[1]}/condor')

import pandas as pd
import config.ana_cff as cfg
from modules.AnaDict import AnaDict

parser = argparse.ArgumentParser(description='Process condor job results and transfer to main file system')
parser.add_argument('-j', dest='jdir', type=str, required=True, help='parse job directory')
#
args = parser.parse_args()

class processResults():

    def __init__(self):
        # check which files need to be transfered
        results = sb.check_output(f'ls {args.jdir}', shell=True).decode().strip()
        all_samples = cfg.All_MC+cfg.Data_samples
        for sample in all_samples:
            # look for non tagged aka non-mini-batch output
            s = re.findall(rf'{sample}_201\d_[a-zA-Z0-9]+.pkl', results) # format sample_year_tag_type.pkl
            for f in s:
                year = re.findall(r'201\d', s[0])[0]
                tag  = f.split('_')[-1].split('.')[0]
                print(f"../files/{year}/{'mc_files' if sample in cfg.All_MC else 'data_files'}/{sample}_{tag}.pkl")
                os.system(
                    f"mv {args.jdir}/{f} ../files/{year}/{'mc_files' if sample in cfg.All_MC else 'data_files'}/{sample}_{tag}.pkl"
                    )
            # now look for tagged aka mini-batch output
            s = re.findall(rf'{sample}_201\d_\d+_[a-zA-Z0-9]+.pkl' , results) # format sample_year_tag_type.pkl, under new format this should always trigger
            if len(s) > 0:
                self.combine_mini_batch(sample, s)

    def combine_mini_batch(self, sample, s_to_comb): 
        # loop of possibilities and join by case
        allowed_types = ['ak4','ak8','rtc','gen','val']
        allowed_years = cfg.Years
        for t in allowed_types:
            for y in allowed_years:
                out_dir = f"../files/{y}/{'mc_files' if sample in cfg.All_MC else 'data_files'}/"
                files = re.findall(rf'{sample}_{y}_\w+_{t}.pkl', ' '.join(s_to_comb))
                if len(files) > 1 :
                    if t == 'val': # vals is in pandas dataframe format
                        df = pd.concat([pd.read_pickle(f'{args.jdir}/{f}') for f in files], ignore_index=True)
                        df.to_pickle(f"{out_dir}{sample}_{t}.pkl")
                    else: # to handle jagged arrays (not sure if this works)
                        a_dict = AnaDict({})
                        [a_dict.update(AnaDict.read_pickle(f'{args.jdir}/{f}')) for f in files]
                        print(f"{out_dir}{sample}_{t}.pkl")
                        a_dict.to_pickle(f"{out_dir}{sample}_{t}.pkl")
                    #
                    [os.system(f'rm {args.jdir}/{f}') for f in files]
                elif len(files) == 1:
                    print(f"{out_dir}{sample}_{t}.pkl")
                    [os.system(f'mv {args.jdir}/{files[0]} {out_dir}{sample}_{t}.pkl')]
                    

if __name__ == '__main__':

    processResults()
